#!/usr/bin/env bash
# Build and publish a tronko reference database on a Jetstream2 VM, from here.
#
#   ops/release/launch.sh --marker CO1_BF3_BR2 --rcrux-tag CO1_BF3_BR2_v1.1.0
#
# It creates a VM that pulls the release image CI pushed for --ref (every commit
# on main has one), runs it, pushes the publish branch and opens its PR, then
# powers off. This script waits, prints the PR and deletes the VM. Merging the PR
# tags the release.
# Closing the laptop is fine: "launch.sh wait NAME" picks the VM up again.
#
#   launch.sh [--ref main] [--flavor m3.2xl] [--volume GB] [--dry-run] OPTIONS
#   launch.sh list           the release VMs and their states
#   launch.sh logs NAME      a VM's console: the build's output so far
#   launch.sh wait NAME      wait for a VM; print its PR and delete it, or show why it failed
#   launch.sh delete NAME    delete a VM and its disk
#
# OPTIONS are scripts/release_tronko.sh's (--set KEY=VALUE, --sp-threshold,
# --bump, ...) and are checked here before anything is created. --config,
# --steps, --work and --tronko-checkout do not apply to a VM.
#
# Needs gcloud, openstack, gh and jq; the js2-worker service account key at
# ~/.rcrux/gcs-key.json (it reads the Jetstream credentials); and
# AWS_ACCESS_KEY_ID / AWS_SECRET_ACCESS_KEY for the Jetstream S3 remotes.
# The VM gets those keys and your GitHub token (GH_TOKEN, else gh's) in its
# user data, and is deleted when it succeeds.
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

wait_for() {
    openstack server show "$1" > /dev/null  # no such VM stops here; a later lost connection only retries
    echo "waiting for $1 (Ctrl-C is safe; '$0 wait $1' picks it up again)"
    until [ "$(openstack server show "$1" -f value -c status)" = SHUTOFF ]; do sleep 60; done
    local result
    result="$(openstack console log show "$1" | tr -d '\r' | grep -o 'TRONKO_RELEASE .*' | tail -1 || true)"
    case "$result" in
        "TRONKO_RELEASE ok "*)
            echo "${result#TRONKO_RELEASE ok }"
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

REF=main DRY_RUN=0 MARKER="" ARGS=()
while [ $# -gt 0 ]; do
    case "$1" in
        --ref) REF="$2"; shift 2 ;;
        --flavor) FLAVOR="$2"; shift 2 ;;
        --volume) VOLUME="$2"; shift 2 ;;
        --dry-run) DRY_RUN=1; shift ;;
        --marker) MARKER="$2"; ARGS+=("$1" "$2"); shift 2 ;;
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
        --config | --steps | --work | --tronko-checkout | --print-config)
            die "$1 does not work here; settings go in --set KEY=VALUE"
            ;;
        -h | --help) usage ;;
        *) ARGS+=("$1"); shift ;;
    esac
done
[ -n "$MARKER" ] || die "--marker is required"
case " ${ARGS[*]} " in *" --rcrux-tag "*) ;; *) die "--rcrux-tag is required" ;; esac
"$HERE/../../scripts/release_tronko.sh" "${ARGS[@]}" --print-config > /dev/null
# Publish refuses a marker main already has without --replace; say so before the build, not after.
case " ${ARGS[*]} " in
    *" --replace "*) ;;
    *) ! gh api "repos/eDNA-Explorer/tronko/contents/databases/$MARKER/tronko.dvc" > /dev/null 2>&1 \
        || die "main already has a tronko release of $MARKER; pass --replace to publish the next version" ;;
esac

# The commit behind a branch, a tag or a short SHA; the image is tagged with it.
SHA="$(gh api "repos/eDNA-Explorer/tronko/commits/$REF" --jq .sha 2> /dev/null)" || die "no tronko commit $REF"

: "${AWS_ACCESS_KEY_ID:?set the Jetstream S3 keys}" "${AWS_SECRET_ACCESS_KEY:?set the Jetstream S3 keys}"
GH_TOKEN="${GH_TOKEN:-$(gh auth token)}"
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
    for key in IMAGE REGISTRY_TOKEN AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY GH_TOKEN GIT_NAME GIT_EMAIL; do
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
