
add_library( angiotk_meshing_volumefromstl
  SHARED
  volumefromstl.cpp
  )
target_link_libraries( angiotk_meshing_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES} )


angiotk_add_component( meshing_centerlines centerlines_main.cpp )
target_link_libraries( meshing_centerlines ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( centerlines.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )

angiotk_add_component( meshing_remeshstl remeshstl_main.cpp )
target_link_libraries( meshing_remeshstl ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
#configure_file( remeshstl.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstl.cfg )

angiotk_add_component( meshing_volumefromstlandcenterlines volumefromstlandcenterlines_main.cpp )
target_link_libraries( meshing_volumefromstlandcenterlines ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
#configure_file( volumefromstlandcenterlines.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstlandcenterlines.cfg )

angiotk_add_component( meshing_volumefromstl volumefromstl_main.cpp )
target_link_libraries( meshing_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( volumefromstl.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstl.cfg )



ExternalData_expand_arguments( angioTkData testArgs DATA{${RemeshBasedOnCenterlines_SOURCE_DIR}/simpleTbifurcation.stl} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${RemeshBasedOnCenterlines_SOURCE_DIR}/simpleTbifurcation_remeshVMTKarea0.005.stl} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${RemeshBasedOnCenterlines_SOURCE_DIR}/simpleTbifurcation.inletoutlet.desc} )
add_dependencies( meshing_volumefromstl angioTkData)
