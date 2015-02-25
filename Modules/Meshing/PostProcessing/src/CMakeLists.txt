
add_library( angiotk_meshing_postprocessing
  SHARED
  postprocessing.cpp
  )
target_link_libraries( angiotk_meshing_postprocessing ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )


angiotk_add_component( meshing_extractsubmesh extractsubmesh_main.cpp )
target_link_libraries( meshing_extractsubmesh ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )
configure_file( extractsubmesh.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_extractsubmesh.cfg )

angiotk_add_component( meshing_meshpartitioner meshpartitioner.cpp )
target_link_libraries( meshing_meshpartitioner ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )
configure_file( meshpartitioner.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_meshpartitioner.cfg )

angiotk_add_component( meshing_postprocessing postprocessing_main.cpp )
target_link_libraries( meshing_postprocessing ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )
configure_file( postprocessing.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_postprocessing.cfg )


ExternalData_expand_arguments( angioTkData testArgs DATA{${PostProcessing_SOURCE_DIR}/simpleTbifurcationFSI.msh} )
add_dependencies(meshing_extractsubmesh angioTkData)
add_dependencies(meshing_meshpartitioner angioTkData)
add_dependencies(meshing_postprocessing angioTkData)
