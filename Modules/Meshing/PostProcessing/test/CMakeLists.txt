
add_test(NAME test_meshing_extractsubmesh COMMAND angiotk_meshing_extractsubmesh
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_extractsubmesh.cfg )

add_test(NAME test_meshing_meshpartitioner COMMAND angiotk_meshing_meshpartitioner
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_meshpartitioner.cfg )

add_test(NAME test_meshing_postprocessing COMMAND angiotk_meshing_pipeline_postprocessing 
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_postprocessing.cfg )
