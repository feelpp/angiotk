
export ANGIOTK_BINARY_DIR=$HOME/AngioTK/angiotk.build/bin/
export MY_CFG_DIR=$HOME/AngioTK/angiotk/Examples/MeshFromMRI/Base2_GEO26/veinous/

#vmtkimagecompose -ifile TOF_vein_cor_GEO26_RORPO_30_1.34_7.nii -i2file TOF_vein_trans_GEO26_RORPO_30_1.34_7.nii -operation max -ofile TOF_vein_fusionmax_GEO26_RORPO_30_1.34_7.mha 
$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage.cfg || exit 1
$ANGIOTK_BINARY_DIR/meshing_centerlines --config-file $MY_CFG_DIR/centerlines.cfg || exit 1
$ANGIOTK_BINARY_DIR/meshing_centerlinesmanager --config-file $MY_CFG_DIR/centerlinesmanager.cfg || exit 1

$ANGIOTK_BINARY_DIR/meshing_imagefromcenterlines --config-file $MY_CFG_DIR/imagefromcenterlines.cfg  || exit 1

$ANGIOTK_BINARY_DIR/meshing_surfacefromimage --config-file $MY_CFG_DIR/surfacefromimage2.cfg || exit 1
$ANGIOTK_BINARY_DIR/meshing_remeshstl --config-file $MY_CFG_DIR/remeshstlgmsh.cfg || exit 1

# $ANGIOTK_BINARY_DIR/meshing_volumefromstlandcenterlines --config-file $MY_CFG_DIR/volumefromstlandcenterlines.cfg || exit 1

