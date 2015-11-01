
export ANGIOTK_BINARY_DIR=$HOME/AngioTK/angiotk.build/bin/
export MY_CFG_DIR=$HOME/AngioTK/angiotk/Examples/MeshFromMRI/Base2_FOA42/

$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage.cfg || exit 1
for CenterlineId in 1 2 3 4 5 6 7 8 9 10
do
    $ANGIOTK_BINARY_DIR/meshing_centerlines --config-file $MY_CFG_DIR/centerlines.cfg \
	--input.pointset.filename=$MY_CFG_DIR/data/pointset$CenterlineId.data \
	--output.directory=angiotk/Base2/FOA42/centerlines/part$CenterlineId || exit 1
done
$ANGIOTK_BINARY_DIR/meshing_centerlinesmanager --config-file $MY_CFG_DIR/centerlinesmanager.cfg || exit 1
$ANGIOTK_BINARY_DIR/meshing_imagefromcenterlines --config-file $MY_CFG_DIR/imagefromcenterlines.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage2.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_remeshstl --config-file $MY_CFG_DIR/remeshstlgmsh.cfg  || exit 1
$ANGIOTK_BINARY_DIR/meshing_volumefromstlandcenterlines --config-file $MY_CFG_DIR/volumefromstlandcenterlines.cfg  || exit 1


