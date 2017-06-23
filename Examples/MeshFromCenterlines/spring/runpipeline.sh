
export ANGIOTK_BINARY_DIR=/usr/local/bin
export MY_CFG_DIR=/feel/AngioTkMeshRessort

$ANGIOTK_BINARY_DIR/angiotk_meshing_centerlinesmanager --config-file $MY_CFG_DIR/centerlinesmanager.cfg || exit 1
$ANGIOTK_BINARY_DIR/angiotk_meshing_imagefromcenterlines --config-file $MY_CFG_DIR/imagefromcenterlines.cfg  || exit 1
$ANGIOTK_BINARY_DIR/angiotk_meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage.cfg || exit 1
$ANGIOTK_BINARY_DIR/angiotk_meshing_remeshstl --config-file $MY_CFG_DIR/remeshstlgmsh.cfg || exit 1
$ANGIOTK_BINARY_DIR/angiotk_meshing_volumefromstlandcenterlines --config-file $MY_CFG_DIR/volumefromstlandcenterlines.cfg  || exit 1
