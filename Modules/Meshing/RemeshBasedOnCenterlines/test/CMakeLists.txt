# Module Meshing
# Label volume to Mesh for CFD

ExternalData_expand_arguments( angioTkData testArgs DATA{${RemeshBasedOnCenterlines_SOURCE_DIR}/simpleTbifurcation.stl} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${RemeshBasedOnCenterlines_SOURCE_DIR}/simpleTbifurcation.inletoutlet.desc} )

add_executable(bloodflowmeshfromSTL bloodflowmeshfromSTL.cpp )
target_link_libraries(bloodflowmeshfromSTL ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})
add_dependencies(bloodflowmeshfromSTL angioTkData)

configure_file( bloodflowmeshfromSTL.cfg ${AngioTk_BINARY_DIR}/bin/bloodflowmeshfromSTL.cfg )
add_test(NAME bloodflowmeshfromSTLTest COMMAND bloodflowmeshfromSTL --config-file=${AngioTk_BINARY_DIR}/bin/bloodflowmeshfromSTL.cfg )
#feelpp_add_application( bloodflowmeshfromSTL SRCS bloodflowmeshfromSTL.cpp INCLUDE_IN_ALL CFG bloodflowmeshfromSTL.cfg )

