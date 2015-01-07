# Module Meshing
# Label volume to Mesh for CFD

add_test(NAME test_meshing_centerlines COMMAND meshing_centerlines
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )

add_test(NAME test_meshing_volumefromstl COMMAND meshing_volumefromstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstl.cfg )
