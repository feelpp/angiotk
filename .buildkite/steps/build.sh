#!/bin/bash

set -eo pipefail
set -x

BRANCH=${BRANCH:-${BUILDKITE_BRANCH:master}}

echo "--- Dockerizing $PROJECT..."
docker run --rm -it -v /var/run/docker.sock:/var/run/docker.sock -e GITHUB_OAUTH  feelpp/feelpp-libs:develop-ubuntu-17.04 \
       sudo \
       PROJECT=${PROJECT} \
       GITHUB_OAUTH=$GITHUB_OAUTH \
       CXX=${CXX} \
       TARGET=${TARGET} \
       FROM=${FROM} \
       BRANCH=${BRANCH} \
       FEELPP_BRANCH=${FEELPP_BRANCH} \
       BUILD_JOBS=${BUILD_JOBS} \
       feelpp_dockerize.sh
