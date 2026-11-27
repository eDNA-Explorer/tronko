# The release VM's boot script. launch.sh puts "#!/bin/bash" and these settings
# above it and passes the result as the server's user data:
#   IMAGE (the release image, tagged with its tronko commit) REGISTRY_TOKEN
#   AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY GIT_NAME GIT_EMAIL
#   RCRUX_COMMIT RCRUX_CURATED_DVC RCRUX_DATA_DVC RCRUX_DVC_CONFIG FORWARD_PRIMER REVERSE_PRIMER
#   (the rcrux tag and its primers, read by launch.sh)
#   ARGS (release_tronko.sh options)
#
# It holds no GitHub credential. tronko is public; the rcrux inputs come from
# the rcrux S3 remote by the hash launch.sh read; the release commit is made
# here and handed back on the console as a git bundle, which launch.sh pushes
# and opens the PR for with your own GitHub login.
#
# Everything it prints goes to the server's console log, which is what
# "launch.sh logs" shows, between "TRONKO_RELEASE start" and "TRONKO_RELEASE ok
# <branch>" or "TRONKO_RELEASE failed"; then the VM powers off, which is what
# launch.sh waits for. TRONKO_MEM lines report the VM's memory use.
set -euo pipefail
exec > >(tee /var/log/tronko-release.log > /dev/console) 2>&1
export HOME=/root
trap 'rc=$?; echo "TRONKO_MEM peak $(cat /run/mem_peak) GiB"; [ "$rc" = 0 ] || echo "TRONKO_RELEASE failed"; poweroff' EXIT
echo "TRONKO_RELEASE start"

# Memory in use (MemTotal - MemAvailable, GiB), sampled every 10 s: a TRONKO_MEM
# line every 10 minutes and the peak on exit, to size the flavor for a marker.
echo 0 > /run/mem_peak
(
    i=0
    while sleep 10; do
        used="$(awk '/^MemTotal:/ {t = $2} /^MemAvailable:/ {a = $2} END {printf "%.1f", (t - a) / 1048576}' /proc/meminfo)"
        awk -v u="$used" -v p="$(cat /run/mem_peak)" 'BEGIN {exit !(u > p)}' && echo "$used" > /run/mem_peak
        (( ++i % 60 )) || echo "TRONKO_MEM $(date -u +%FT%TZ) used $used GiB, peak $(cat /run/mem_peak) GiB"
    done
) &

mkdir -p /work/pointer
git clone --quiet https://github.com/eDNA-Explorer/tronko.git /work/release
printf '%s\n' "$RCRUX_COMMIT" > /work/pointer/commit
printf '%s\n' "$RCRUX_CURATED_DVC" > /work/pointer/curated.dvc
printf '%s\n' "$RCRUX_DATA_DVC" > /work/pointer/data.dvc
printf '%s\n' "$RCRUX_DVC_CONFIG" > /work/pointer/dvc_config
printf '%s %s\n' "$FORWARD_PRIMER" "$REVERSE_PRIMER" > /work/pointer/primers
echo "$REGISTRY_TOKEN" | docker login -u oauth2accesstoken --password-stdin "https://${IMAGE%%/*}"
docker pull --quiet "$IMAGE"

# --ipc=host: tronko-assign (the assignment check) stages its BWA-MEM3 index in
# /dev/shm, and a container's own is 64 MB.
docker run --rm --ipc=host -v /work:/work \
    -e AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" -e AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    -e GIT_AUTHOR_NAME="$GIT_NAME" -e GIT_AUTHOR_EMAIL="$GIT_EMAIL" \
    -e GIT_COMMITTER_NAME="$GIT_NAME" -e GIT_COMMITTER_EMAIL="$GIT_EMAIL" \
    "$IMAGE" --work /work/build --tronko-checkout /work/release \
    --rcrux-pointer /work/pointer "${ARGS[@]}"

# publish left the checkout on its branch, one commit on origin/main: the PR's
# title and, below it, the QC report. A few KB, so the console carries it.
branch="$(git -C /work/release branch --show-current)"
git -C /work/release bundle create --quiet /work/release.bundle "origin/main..$branch"
echo "TRONKO_BUNDLE $branch $(sha256sum /work/release.bundle | cut -d' ' -f1)"
base64 -w 76 /work/release.bundle
echo "TRONKO_BUNDLE end"
echo "TRONKO_RELEASE ok $branch"
