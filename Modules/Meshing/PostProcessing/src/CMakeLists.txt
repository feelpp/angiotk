
add_library( angiotk_meshing_postprocessing
  SHARED
  postprocessing.cpp
  )
target_link_libraries( angiotk_meshing_postprocessing Feelpp::feelpp)
install(TARGETS angiotk_meshing_postprocessing LIBRARY DESTINATION lib/ COMPONENT Libs)

angiotk_add_component( angiotk_meshing_extractsubmesh extractsubmesh_main.cpp )
target_link_libraries( angiotk_meshing_extractsubmesh ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )

angiotk_add_component( angiotk_meshing_meshpartitioner meshpartitioner.cpp )
target_link_libraries( angiotk_meshing_meshpartitioner ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )

angiotk_add_component( angiotk_meshing_pipeline_postprocessing postprocessing_main.cpp )
target_link_libraries( angiotk_meshing_pipeline_postprocessing ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_postprocessing )


