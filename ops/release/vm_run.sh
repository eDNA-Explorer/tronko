# The release VM's boot script. launch.sh puts "#!/bin/bash" and these settings
# above it and passes the result as the server's user data:
#   IMAGE REGISTRY_TOKEN AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY GH_TOKEN
#   GIT_NAME GIT_EMAIL ARGS (release_tronko.sh options)
#
# Everything it prints goes to the server's console log, which is what
# "launch.sh logs" shows. The last line is "TRONKO_RELEASE ok <PR url>" or
# "TRONKO_RELEASE failed"; then the VM powers off, which is what launch.sh
# waits for.
set -euo pipefail
exec > >(tee /var/log/tronko-release.log > /dev/console) 2>&1
export HOME=/root
trap '[ "$?" = 0 ] || echo "TRONKO_RELEASE failed"; poweroff' EXIT

GITHUB_REPO=eDNA-Explorer/tronko

echo "$REGISTRY_TOKEN" | docker login -u oauth2accesstoken --password-stdin "https://${IMAGE%%/*}"
docker pull --quiet "$IMAGE"
mkdir -p /work
git clone --quiet "https://x-access-token:$GH_TOKEN@github.com/$GITHUB_REPO.git" /work/release

docker run --rm -v /work:/work \
    -e AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID" -e AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY" \
    -e GIT_AUTHOR_NAME="$GIT_NAME" -e GIT_AUTHOR_EMAIL="$GIT_EMAIL" \
    -e GIT_COMMITTER_NAME="$GIT_NAME" -e GIT_COMMITTER_EMAIL="$GIT_EMAIL" \
    "$IMAGE" --work /work/build --tronko-checkout /work/release --push "${ARGS[@]}"

# publish left the checkout on the branch it pushed; its one commit is the title.
branch="$(git -C /work/release branch --show-current)"
title="$(git -C /work/release log -1 --format=%s)"
pr="$(curl -fsS https://api.github.com/repos/$GITHUB_REPO/pulls \
    -H "Authorization: Bearer $GH_TOKEN" -H "Accept: application/vnd.github+json" \
    -d "{\"title\": \"$title\", \"head\": \"$branch\", \"base\": \"main\"}" \
    | sed -n 's|.*"html_url": "\(.*/pull/[0-9]*\)".*|\1|p' | head -1)"
echo "TRONKO_RELEASE ok $pr"
