
export ANGIOTK_BINARY_DIR=$HOME/AngioTK/angiotk.build/bin/
export MY_CFG_DIR=$HOME/AngioTK/angiotk/Examples/MeshFromMRI/Base1_000/

$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage.cfg  || exit 1

$ANGIOTK_BINARY_DIR/meshing_centerlines --config-file $MY_CFG_DIR/centerlines.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_centerlines --config-file $MY_CFG_DIR/centerlines.cfg \
    --input.pointpair.filename=\$cfgdir/data/pointpair2.data --output.directory=angiotk/Base1/000/centerlines/part2  || exit 1
$ANGIOTK_BINARY_DIR/meshing_centerlines --config-file $MY_CFG_DIR/centerlines.cfg \
    --input.geo-centerlines.filename=\$cfgdir/data/geocenterlines.data --output.directory=angiotk/Base1/000/centerlines/part3 \
    --post-process.convert-centerlines=1 --convert-centerlines.threshold-radius.min=1.5  || exit 1

# $ANGIOTK_BINARY_DIR/meshing_centerlinesmanager --config-file $MY_CFG_DIR/centerlinesfusion.cfg
$ANGIOTK_BINARY_DIR/meshing_centerlinesmanager --config-file $MY_CFG_DIR/centerlinesmanager.cfg  || exit 1

$ANGIOTK_BINARY_DIR/meshing_imagefromcenterlines --config-file $MY_CFG_DIR/imagefromcenterlines.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage2.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_remeshstl --config-file $MY_CFG_DIR/remeshstlgmsh.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_volumefromstlandcenterlines --config-file $MY_CFG_DIR/volumefromstlandcenterlines.cfg || exit 1

