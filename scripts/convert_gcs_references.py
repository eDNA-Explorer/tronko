#!/usr/bin/env python3
"""
Convert reference_tree.txt.gz files from GCS to zstd-compressed .trkb format.

This script:
1. Lists all marker folders in gs://edna-project-files-staging/CruxV2/
2. Downloads reference_tree.txt.gz for each marker
3. Converts to zstd-compressed .trkb using tronko-convert
4. Uploads the .trkb file back to GCS
5. Cleans up local temporary files

Usage:
    python convert_gcs_references.py [--dry-run] [--marker MARKER] [--compression-level LEVEL]

Requirements:
    pip install google-cloud-storage
"""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

from google.cloud import storage


BUCKET_NAME = "edna-project-files-staging"
BASE_PREFIX = "CruxV2"
TRONKO_DATE = "2023-04-07"
TRONKO_CONVERT_PATH = Path(__file__).parent.parent / "tronko-convert" / "tronko-convert"


def list_markers(client: storage.Client, bucket_name: str, base_prefix: str) -> list[str]:
    """List all marker folders under the base prefix."""
    bucket = client.bucket(bucket_name)

    # List with delimiter to get "folders"
    blobs = bucket.list_blobs(prefix=f"{base_prefix}/", delimiter="/")

    # Consume the iterator to populate prefixes
    list(blobs)

    markers = []
    for prefix in blobs.prefixes:
        # prefix looks like "CruxV2/16S_Bacteria/"
        marker = prefix.rstrip("/").split("/")[-1]
        markers.append(marker)

    return sorted(markers)


def download_file(client: storage.Client, bucket_name: str, blob_path: str, local_path: Path) -> bool:
    """Download a file from GCS. Returns True if successful."""
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)

    if not blob.exists():
        return False

    blob.download_to_filename(str(local_path))
    return True


def upload_file(client: storage.Client, bucket_name: str, local_path: Path, blob_path: str) -> None:
    """Upload a file to GCS."""
    bucket = client.bucket(bucket_name)
    blob = bucket.blob(blob_path)
    blob.upload_from_filename(str(local_path))


def convert_reference(input_path: Path, output_path: Path, compression_level: int = 19) -> bool:
    """Convert reference_tree.txt.gz to zstd-compressed .trkb using tronko-convert."""
    cmd = [
        str(TRONKO_CONVERT_PATH),
        "-i", str(input_path),
        "-o", str(output_path),
        "-c", str(compression_level),
        "-v"
    ]

    result = subprocess.run(cmd, capture_output=True, text=True)

    if result.returncode != 0:
        print(f"  ERROR: tronko-convert failed:")
        print(f"  stdout: {result.stdout}")
        print(f"  stderr: {result.stderr}")
        return False

    return True


def process_marker(
    client: storage.Client,
    marker: str,
    compression_level: int = 19,
    dry_run: bool = False
) -> tuple[bool, str]:
    """
    Process a single marker: download, convert, upload, cleanup.

    Returns (success, message).
    """
    source_path = f"{BASE_PREFIX}/{marker}/{TRONKO_DATE}/tronko/reference_tree.txt.gz"
    dest_path = f"{BASE_PREFIX}/{marker}/{TRONKO_DATE}/tronko/reference_tree.trkb"

    if dry_run:
        return True, f"Would convert gs://{BUCKET_NAME}/{source_path} -> gs://{BUCKET_NAME}/{dest_path}"

    with tempfile.TemporaryDirectory() as tmpdir:
        tmpdir = Path(tmpdir)
        local_source = tmpdir / "reference_tree.txt.gz"
        local_dest = tmpdir / "reference_tree.trkb"

        # Download
        print(f"  Downloading {source_path}...")
        if not download_file(client, BUCKET_NAME, source_path, local_source):
            return False, f"Source file not found: gs://{BUCKET_NAME}/{source_path}"

        source_size = local_source.stat().st_size / (1024 * 1024)
        print(f"  Downloaded {source_size:.1f} MB")

        # Convert
        print(f"  Converting to zstd-compressed .trkb (level {compression_level})...")
        if not convert_reference(local_source, local_dest, compression_level):
            return False, "Conversion failed"

        dest_size = local_dest.stat().st_size / (1024 * 1024)
        ratio = source_size / dest_size if dest_size > 0 else 0
        print(f"  Converted: {dest_size:.1f} MB ({ratio:.1f}x compression vs gzipped text)")

        # Upload
        print(f"  Uploading to {dest_path}...")
        upload_file(client, BUCKET_NAME, local_dest, dest_path)
        print(f"  Upload complete")

        # Cleanup happens automatically when tmpdir context exits

    return True, f"Successfully converted {marker}: {source_size:.1f} MB -> {dest_size:.1f} MB ({ratio:.1f}x)"


def main():
    parser = argparse.ArgumentParser(
        description="Convert GCS reference_tree.txt.gz files to zstd-compressed .trkb"
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        help="List markers and show what would be done without actually converting"
    )
    parser.add_argument(
        "--marker",
        type=str,
        help="Process only this specific marker (default: all markers)"
    )
    parser.add_argument(
        "--compression-level",
        type=int,
        default=19,
        choices=range(1, 20),
        metavar="1-19",
        help="Zstd compression level (default: 19, max compression)"
    )
    parser.add_argument(
        "--list-only",
        action="store_true",
        help="Only list available markers, don't process anything"
    )

    args = parser.parse_args()

    # Verify tronko-convert exists
    if not TRONKO_CONVERT_PATH.exists():
        print(f"ERROR: tronko-convert not found at {TRONKO_CONVERT_PATH}")
        print("Please build tronko-convert first: cd tronko-convert && make")
        sys.exit(1)

    # Initialize GCS client
    client = storage.Client()

    # List markers
    print(f"Listing markers in gs://{BUCKET_NAME}/{BASE_PREFIX}/...")
    markers = list_markers(client, BUCKET_NAME, BASE_PREFIX)
    print(f"Found {len(markers)} markers: {', '.join(markers)}")
    print()

    if args.list_only:
        return

    # Filter to specific marker if requested
    if args.marker:
        if args.marker not in markers:
            print(f"ERROR: Marker '{args.marker}' not found in available markers")
            sys.exit(1)
        markers = [args.marker]

    # Process markers
    results = []
    for i, marker in enumerate(markers, 1):
        print(f"[{i}/{len(markers)}] Processing {marker}...")
        success, message = process_marker(
            client,
            marker,
            compression_level=args.compression_level,
            dry_run=args.dry_run
        )
        results.append((marker, success, message))
        print(f"  {message}")
        print()

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)

    succeeded = [r for r in results if r[1]]
    failed = [r for r in results if not r[1]]

    print(f"Succeeded: {len(succeeded)}/{len(results)}")
    if failed:
        print(f"Failed: {len(failed)}")
        for marker, _, message in failed:
            print(f"  - {marker}: {message}")

    sys.exit(0 if not failed else 1)


if __name__ == "__main__":
    main()
