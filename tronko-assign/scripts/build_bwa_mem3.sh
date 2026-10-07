#!/usr/bin/env bash
# Build the pinned BWA-MEM3 that tronko-assign runs as its aligner (a separate program; it is
# never linked into tronko-assign). Used by `make bwa-mem3` in tronko-assign/.
#
#   build_bwa_mem3.sh [build-dir]       default build dir: tronko-assign/build/bwa-mem3
#
# Clones https://github.com/fg-labs/bwa-mem3 at tag v0.14.0 with its submodules into
# <build-dir>/src (a fresh clone when the directory is absent), refuses to build unless the
# checkout is commit 5c1d5e39e0b05fe3cf960ed0dc6086af59b3e870, builds it unchanged (no patch to
# its sources) and copies the binary to <build-dir>/bwa-mem3.
#
# Compiler: BWA-MEM3 requires clang 19 or GCC 15. The default is CC=clang-19 CXX=clang++-19;
# set CC/CXX to use GCC 15 instead. The shipped binary must be built with clang 19 or GCC 15.
# For development only, ALLOW_UNSUPPORTED_COMPILER=1 lets BWA-MEM3's Makefile accept an older
# compiler (for example the gcc 13 of Ubuntu 24.04): CC=gcc CXX=g++ ALLOW_UNSUPPORTED_COMPILER=1.
#
# Environment: BWA_MEM3_JOBS (parallel make jobs, default nproc), CC, CXX,
# ALLOW_UNSUPPORTED_COMPILER, BWA_MEM3_MAKE_ARGS (extra arguments to BWA-MEM3's make, e.g. a
# target tier; empty by default).
set -euo pipefail

BWA_MEM3_REPO=https://github.com/fg-labs/bwa-mem3
BWA_MEM3_TAG=v0.14.0
BWA_MEM3_COMMIT=5c1d5e39e0b05fe3cf960ed0dc6086af59b3e870

HERE=$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)
BUILD=${1:-"$HERE/build/bwa-mem3"}
SRC=$BUILD/src
CC=${CC:-clang-19}
CXX=${CXX:-clang++-19}
JOBS=${BWA_MEM3_JOBS:-$(nproc)}

mkdir -p "$BUILD"
if [[ ! -d $SRC/.git ]]; then
	rm -rf "$SRC"
	git clone --quiet --branch "$BWA_MEM3_TAG" --recurse-submodules "$BWA_MEM3_REPO" "$SRC"
fi
got=$(git -C "$SRC" rev-parse HEAD)
if [[ $got != "$BWA_MEM3_COMMIT" ]]; then
	echo "build_bwa_mem3.sh: $SRC is at $got, expected $BWA_MEM3_TAG = $BWA_MEM3_COMMIT" >&2
	exit 1
fi
git -C "$SRC" submodule update --init --recursive --quiet

command -v "$CC" >/dev/null || { echo "build_bwa_mem3.sh: compiler $CC not found (set CC/CXX; see the header)" >&2; exit 1; }
command -v "$CXX" >/dev/null || { echo "build_bwa_mem3.sh: compiler $CXX not found (set CC/CXX; see the header)" >&2; exit 1; }

extra=()
[[ ${ALLOW_UNSUPPORTED_COMPILER:-0} == 1 ]] && extra+=(ALLOW_UNSUPPORTED_COMPILER=1)
# shellcheck disable=SC2206
extra+=(${BWA_MEM3_MAKE_ARGS:-})
make -C "$SRC" -j"$JOBS" CC="$CC" CXX="$CXX" "${extra[@]}"
cp "$SRC/bwa-mem3" "$BUILD/bwa-mem3"
{
	echo "repo $BWA_MEM3_REPO"
	echo "tag $BWA_MEM3_TAG"
	echo "commit $got"
	echo "cc $("$CC" --version | head -1)"
	echo "cxx $("$CXX" --version | head -1)"
	echo "allow_unsupported_compiler ${ALLOW_UNSUPPORTED_COMPILER:-0}"
	echo "sha256 $(sha256sum "$BUILD/bwa-mem3" | cut -d' ' -f1)"
} >"$BUILD/build-info.txt"
cat "$BUILD/build-info.txt"
