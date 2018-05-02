#!/bin/bash

set -eo pipefail

run="$(basename "$0")"
ANGIOTK_DIR=/usr/local/bin
CFG_DIR=$PWD

usage() {
    echo >&2 "usage: $run "
    echo >&2 "              [--angiotk angiotk directory, default: $DEFAULT_ANGIOTK_DIR]"
    echo >&2 "              [--cfg configuration files for angiotk, default: $CFG_DIR]"
    echo >&2 "   ie: $run --angiotk /usr/local/bin --cfg \$PWD"
    exit 1
}

while [ -n "$1" ]; do
    case "$1" in
        --angiotk) ANGIOTK_DIR="$2" ; shift 2 ;;
        --cfg) CFG_DIR="$2" ; shift 2 ;;
        -h|--help) usage ;;
        --) shift ; break ;;
    esac
done
echo " . angiotk tools: $ANGIOTK_DIR"
echo " . Configuration files for angiotk pipeline: $CFG_DIR"
$ANGIOTK_DIR/angiotk_meshing_centerlinesmanager --config-file $CFG_DIR/centerlinesmanager.cfg || exit 1
$ANGIOTK_DIR/angiotk_meshing_imagefromcenterlines --config-file $CFG_DIR/imagefromcenterlines.cfg  || exit 1
$ANGIOTK_DIR/angiotk_meshing_surfacefromimage --config-file $CFG_DIR/surfacefromimage.cfg || exit 1
$ANGIOTK_DIR/angiotk_meshing_remeshstl --config-file $CFG_DIR/remeshstlgmsh.cfg || exit 1
$ANGIOTK_DIR/angiotk_meshing_volumefromstlandcenterlines --config-file $CFG_DIR/volumefromstlandcenterlines.cfg  || exit 1
