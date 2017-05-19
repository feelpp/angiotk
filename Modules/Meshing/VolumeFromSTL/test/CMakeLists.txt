# Module Meshing
# Label volume to Mesh for CFD

add_test(NAME test_meshing_centerlines COMMAND angiotk_meshing_centerlines
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )

add_test(NAME test_meshing_remeshstlvmtk COMMAND angiotk_meshing_remeshstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlvmtk.cfg )

add_test(NAME test_meshing_remeshstlgmsh COMMAND angiotk_meshing_remeshstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlgmsh.cfg )

add_test(NAME test_meshing_volumefromstlandcenterlines COMMAND angiotk_meshing_volumefromstlandcenterlines
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstlandcenterlines.cfg )

add_test(NAME test_meshing_volumefromstl COMMAND angiotk_meshing_pipeline_volumefromstl
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstl.cfg )
