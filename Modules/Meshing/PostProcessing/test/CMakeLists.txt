
add_test(NAME test_meshing_extractsubmesh COMMAND meshing_extractsubmesh
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_extractsubmesh.cfg )

add_test(NAME test_meshing_meshpartitioner COMMAND meshing_meshpartitioner
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_meshpartitioner.cfg )

add_test(NAME test_meshing_postprocessing COMMAND meshing_postprocessing 
  --config-file=${AngioTk_BINARY_DIR}/cfg/meshing/test_postprocessing.cfg )
