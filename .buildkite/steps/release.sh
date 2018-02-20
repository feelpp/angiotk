#! /bin/bash

set -eo pipefail
set -x

if [ -v DOCKER_PASSWORD -a -v DOCKER_LOGIN ]; then
    docker login --username="${DOCKER_LOGIN}" --password="${DOCKER_PASSWORD}";
fi    

BRANCH=${BRANCH:-${BUILDKITE_BRANCH:master}}

mkdir -p tools/scripts/buildkite/
mkdir -p cmake/modules/

docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs \
       cat /usr/local/share/feelpp/scripts/release.sh | dos2unix > tools/scripts/buildkite/release.sh
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs \
       cat /usr/local/share/feelpp/scripts/list.sh | dos2unix > tools/scripts/buildkite/list.sh
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs \
       cat /usr/local/share/feelpp/scripts/common.sh | dos2unix > tools/scripts/buildkite/common.sh
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock feelpp/feelpp-libs \
       cat /usr/local/share//feelpp/feel/cmake/modules/feelpp.version.cmake | dos2unix > cmake/modules/feelpp.version.cmake

chmod u+x tools/scripts/buildkite/release.sh
chmod u+x tools/scripts/buildkite/list.sh
chmod u+x tools/scripts/buildkite/common.sh

tools/scripts/buildkite/release.sh -t ${TARGET} -b ${BRANCH} -- ${PROJECT}


