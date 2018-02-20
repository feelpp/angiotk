#!/bin/bash

set -eo pipefail
set -x

BRANCH=${BRANCH:-${BUILDKITE_BRANCH:master}}

echo "--- Dockerizing $PROJECT..."
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock -e GITHUB_OAUTH  feelpp/feelpp-libs \
       sudo \
       GITHUB_OAUTH=$GITHUB_OAUTH \
       CXX=${CXX} \
       feelpp_dockerize.sh -p=${PROJECT} -t=${TARGET} -c="${CMAKE_FLAGS}" -f=${FROM} -j=${BUILD_JOBS} -b=${BRANCH} --feelpp-branch=${FEELPP_BRANCH}
