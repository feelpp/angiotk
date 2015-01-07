# Module Meshing
# Label volume to Mesh for CFD

add_test(NAME test_meshing_centerlines COMMAND meshing_centerlines
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )

add_test(NAME test_meshing_remeshstlvmtk COMMAND meshing_remeshstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlvmtk.cfg )

add_test(NAME test_meshing_remeshstlgmsh COMMAND meshing_remeshstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlgmsh.cfg )

add_test(NAME test_meshing_volumefromstlandcenterlines COMMAND meshing_volumefromstlandcenterlines
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstlandcenterlines.cfg )

add_test(NAME test_meshing_volumefromstl COMMAND meshing_volumefromstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstl.cfg )
