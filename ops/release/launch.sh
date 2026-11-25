#!/usr/bin/env bash
# Build and publish a tronko reference database on a Jetstream2 VM, from here.
#
#   ops/release/launch.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.0.0
#
# It creates a VM that pulls the release image CI pushed for --ref (every commit
# on main has one), runs it, pushes the database to the DVC remote and commits
# the release, then powers off. This script waits, pushes the publish branch and
# opens its PR with your git and gh, and deletes the VM. Merging the PR tags the
# release. Closing the laptop is fine: "launch.sh wait NAME" picks the VM up again.
#
#   launch.sh [--ref main] [--flavor m3.2xl] [--volume GB] [--dry-run]
#             [--forward-primer SEQ --reverse-primer SEQ] OPTIONS
#   launch.sh list           the release VMs and their states
#   launch.sh logs NAME      a VM's console: the build's output so far
#   launch.sh wait NAME      wait for a VM; open its PR and delete it, or show why it failed
#   launch.sh delete NAME    delete a VM and its disk
#
# OPTIONS are scripts/release_tronko.sh's (--set KEY=VALUE, --sp-threshold,
# --bump, ...) and are checked here before anything is created. The primers
# (for the assignment check) are read from the rcrux release; the two primer
# options give them for a release that records none, or replace them. --config,
# --steps, --work and --tronko-checkout do not apply to a VM.
#
# Needs gcloud, openstack, gh, git and jq; the js2-worker service account key at
# ~/.rcrux/gcs-key.json (it reads the Jetstream credentials); and
# AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY for the Jetstream S3 remotes.
# The VM gets those keys in its user data and no GitHub credential: this script
# reads the rcrux tag for it, and pushes and opens the PR itself. The VM is
# deleted when it succeeds.
set -euo pipefail

usage() { awk 'NR > 1 && /^#/ { sub(/^# ?/, ""); print; next } NR > 1 { exit }' "$0"; exit "${1:-0}"; }
die() { echo "error: $*" >&2; exit 1; }

HERE="$(cd "$(dirname "$0")" && pwd)"
FLAVOR="${FLAVOR:-m3.2xl}"
VOLUME="${VOLUME:-300}"
VM_IMAGE="${VM_IMAGE:-Featured-Ubuntu22}"  # has Docker
NETWORK="${NETWORK:-auto_allocated_network}"
REGISTRY=northamerica-northeast2-docker.pkg.dev/edna-explorer-shared/jetstream2-images

delete() {
    local volumes
    volumes="$(openstack server show "$1" -f json -c volumes_attached | jq -r '.volumes_attached[].id')"
    openstack server delete --wait "$1"
    [ -z "$volumes" ] || openstack volume delete $volumes
}

# Push the release commit the VM printed as a bundle, and open (or find) its PR.
publish() {
    local console="$1" branch sha bundle repo
    read -r branch sha < <(grep '^TRONKO_BUNDLE [^ ]* [0-9a-f]\{64\}$' <<< "$console" | tail -1 | cut -d' ' -f2-) || true
    [ -n "${branch:-}" ] || die "no release bundle on the console"
    bundle="$(mktemp)"
    awk '/^TRONKO_BUNDLE [^ ]* [0-9a-f]+$/ { body = 1; next } /^TRONKO_BUNDLE end$/ { body = 0 }
        body && /^[A-Za-z0-9+\/=]+$/' <<< "$console" | base64 --decode > "$bundle"
    [ "$( { sha256sum 2> /dev/null || shasum -a 256; } < "$bundle" | cut -d' ' -f1)" = "$sha" ] \
        || die "the release bundle came back damaged; '$0 wait NAME' reads the console again"
    repo="$HERE/../.."
    git -C "$repo" fetch --quiet origin main  # the bundle's base
    git -C "$repo" fetch --quiet "$bundle" "refs/heads/$branch"
    git -C "$repo" push --quiet origin "FETCH_HEAD:refs/heads/$branch"
    gh pr view "$branch" --repo eDNA-Explorer/tronko --json url --jq .url 2> /dev/null \
        || gh pr create --repo eDNA-Explorer/tronko --base main --head "$branch" \
            --title "$(git -C "$repo" log -1 --format=%s FETCH_HEAD)" \
            --body "$(git -C "$repo" log -1 --format=%b FETCH_HEAD)"
    rm -f "$bundle"
}

wait_for() {
    openstack server show "$1" > /dev/null  # no such VM stops here; a later lost connection only retries
    echo "waiting for $1 (Ctrl-C is safe; '$0 wait $1' picks it up again)"
    until [ "$(openstack server show "$1" -f value -c status)" = SHUTOFF ]; do sleep 60; done
    local console result
    console="$(openstack console log show "$1" | tr -d '\r')"
    result="$(grep -o 'TRONKO_RELEASE .*' <<< "$console" | tail -1 || true)"
    case "$result" in
        "TRONKO_RELEASE ok "*)
            publish "$console"
            delete "$1"
            ;;
        *)
            # The end of the boot script's output, without the shutdown that follows it.
            # (OpenStack keeps only the console's last 100 KiB, so its start may be gone.)
            openstack console log show "$1" | tr -d '\r' | sed '/TRONKO_RELEASE failed/q' | tail -40 || true
            die "$1 failed and is kept, powered off; '$0 delete $1' removes it"
            ;;
    esac
}

case "${1:-}" in -h | --help | "") usage ;; esac

export CLOUDSDK_AUTH_CREDENTIAL_FILE_OVERRIDE="${GCP_KEY:-$HOME/.rcrux/gcs-key.json}"
OPENSTACK_CREDENTIALS="$(gcloud secrets versions access latest --secret=openstack-credentials --project=edna-explorer-shared)"
eval "$OPENSTACK_CREDENTIALS"

case "$1" in
    list) openstack server list --name '^tronko-release-' -c Name -c Status; exit ;;
    logs) openstack console log show "${2:?NAME}"; exit ;;
    wait) wait_for "${2:?NAME}"; exit ;;
    delete) delete "${2:?NAME}"; exit ;;
esac

REF=main DRY_RUN=0 MARKER="" RCRUX_TAG="" FORWARD_PRIMER="" REVERSE_PRIMER="" ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --ref) REF="$2"; shift 2 ;;
        --flavor) FLAVOR="$2"; shift 2 ;;
        --volume) VOLUME="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --marker) MARKER="$2"; ARGS+=("$1" "$2"); shift 2 ;;
        --rcrux-tag) RCRUX_TAG="$2"; ARGS+=("$1" "$2"); shift 2 ;;
        --forward-primer) FORWARD_PRIMER="$2"; shift 2 ;;
        --reverse-primer) REVERSE_PRIMER="$2"; shift 2 ;;
        # release_tronko.sh checks these two only when it publishes, after the build.
        --bump)
            case "$2" in major | minor | patch) ;; *) die "--bump must be major, minor or patch" ;; esac
            ARGS+=("$1" "$2"); shift 2
            ;;
        --tronko-version)
            [[ "$2" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || die "--tronko-version must look like v1.0.0"
            ARGS+=("$1" "$2"); shift 2
            ;;
        # The VM cannot read a local file, always runs every step, and sets its own paths.
        --config | --steps | --work | --tronko-checkout | --print-config | --rcrux-repo | --rcrux-pointer | --push)
            die "$1 does not work here; settings go in --set KEY=VALUE"
            ;;
        -h | --help) usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ -n "$MARKER" ] || die "--marker is required"
[ -n "$RCRUX_TAG" ] || die "--rcrux-tag is required"
"$HERE/../../scripts/release_tronko.sh" "${ARGS[@]}" --print-config > /dev/null
# Publish refuses a marker main already has without --replace; say so before the build, not after.
case " ${ARGS[*]} " in
    *" --replace "*) ;;
    *) ! gh api "repos/eDNA-Explorer/tronko/contents/databases/$MARKER/tronko.dvc" > /dev/null 2>&1 \
        || die "main already has a tronko release of $MARKER; pass --replace to publish the next version" ;;
esac

# The rcrux tag, read here with your GitHub access so the VM needs none: its
# commit, the marker's curated.dvc (the hash of the inputs) and the remote's settings.
RCRUX_COMMIT="$(gh api "repos/eDNA-Explorer/rcrux-py/commits/$RCRUX_TAG" --jq .sha 2> /dev/null)" \
    || die "no rcrux-py tag $RCRUX_TAG"
rcrux_file() { gh api "repos/eDNA-Explorer/rcrux-py/contents/$1?ref=$RCRUX_COMMIT" -H "Accept: application/vnd.github.raw"; }
RCRUX_CURATED_DVC="$(rcrux_file "databases/$MARKER/curated.dvc")" || die "$RCRUX_TAG has no databases/$MARKER/curated.dvc"
RCRUX_DVC_CONFIG="$(rcrux_file .dvc/config)"
if [ -z "$FORWARD_PRIMER" ] || [ -z "$REVERSE_PRIMER" ]; then
    recorded="$(python3 "$HERE/../../scripts/rcrux_primers.py" "$MARKER" \
        <(rcrux_file "databases/$MARKER/build.log" 2> /dev/null || true) \
        <(rcrux_file src/rcrux/config.py 2> /dev/null || true))" \
        || die "$RCRUX_TAG records no primers; give --forward-primer and --reverse-primer"
    FORWARD_PRIMER="${FORWARD_PRIMER:-${recorded% *}}" REVERSE_PRIMER="${REVERSE_PRIMER:-${recorded#* }}"
fi

# The commit behind a branch, a tag or a short SHA; the image is tagged with it.
SHA="$(gh api "repos/eDNA-Explorer/tronko/commits/$REF" --jq .sha 2> /dev/null)" || die "no tronko commit $REF"

: "${AWS_ACCESS_KEY_ID:?set the Jetstream S3 keys}" "${AWS_SECRET_ACCESS_KEY:?set the Jetstream S3 keys}"
IMAGE="$REGISTRY/tronko-release:$SHA"
REGISTRY_TOKEN="$(gcloud auth print-access-token)"  # js2-worker's, to pull the image; lasts an hour
curl -fsS -o /dev/null -H "Authorization: Bearer $REGISTRY_TOKEN" \
    "https://${REGISTRY%%/*}/v2/${REGISTRY#*/}/tronko-release/manifests/$SHA" 2> /dev/null \
    || die "no image $IMAGE; build it: gh workflow run build-push-tronko-release.yml -f ref=$SHA"
# The release commit is made as you.
GIT_NAME="$(git config user.name)" && GIT_EMAIL="$(git config user.email)" || die "set git user.name and user.email"

USER_DATA="$(mktemp)"
trap 'rm -f "$USER_DATA"' EXIT
{
    echo '#!/bin/bash'
    for key in IMAGE REGISTRY_TOKEN AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY GIT_NAME GIT_EMAIL \
            RCRUX_COMMIT RCRUX_CURATED_DVC RCRUX_DVC_CONFIG FORWARD_PRIMER REVERSE_PRIMER; do
        printf '%s=%q\n' "$key" "${!key}"
    done
    printf 'ARGS=(%s)\n' "$(printf '%q ' "${ARGS[@]}")"
    cat "$HERE/vm_run.sh"
} > "$USER_DATA"

NAME="tronko-release-$(echo "$MARKER" | tr 'A-Z_' 'a-z-')-$(date +%m%d%H%M)"
CREATE=(openstack server create "$NAME" --flavor "$FLAVOR" --image "$VM_IMAGE" --network "$NETWORK"
    --boot-from-volume "$VOLUME" --user-data "$USER_DATA" --wait)
[ -z "${KEY_NAME:-}" ] || CREATE+=(--key-name "$KEY_NAME")  # only to ssh in from inside Jetstream
if [ "$DRY_RUN" = 1 ]; then
    echo "would run: ${CREATE[*]}"
    exit
fi
"${CREATE[@]}" > /dev/null
wait_for "$NAME"
