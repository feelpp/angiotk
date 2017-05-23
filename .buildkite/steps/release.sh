#! /bin/bash

set -eo pipefail
set -x

if [ -v DOCKER_PASSWORD -a -v DOCKER_LOGIN ]; then
    docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}";
fi    

BRANCH=${BRANCH:-${BUILDKITE_BRANCH:master}}

mkdir -p tools/scripts/buildkite/

docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs:develop-ubuntu-17.04 \
       cat /usr/local/share/feelpp/scripts/release.sh | dos2unix > tools/scripts/buildkite/release.sh
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs:develop-ubuntu-17.04 \
       cat /usr/local/share/feelpp/scripts/list.sh | dos2unix > tools/scripts/buildkite/list.sh

chmod u+x tools/scripts/buildkite/release.sh
chmod u+x tools/scripts/buildkite/list.sh

tools/scripts/buildkite/release.sh -t ${TARGET} -b ${BRANCH} -- ${PROJECT}


