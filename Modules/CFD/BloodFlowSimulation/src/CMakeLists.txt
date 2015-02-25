
angiotk_add_component( cfd_bloodflow bloodflow_main.cpp )
target_link_libraries( cfd_bloodflow ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})
#configure_file( centerlines.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )
