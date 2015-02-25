
add_library( angiotk_meshing_volumefromstl
  SHARED
  volumefromstl.cpp
  AngioTkCenterlineField.cpp
  )
target_link_libraries( angiotk_meshing_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} )


angiotk_add_component( meshing_createinletoutletdesc createinletoutletdesc_main.cpp )
target_link_libraries( meshing_createinletoutletdesc ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_centerlines centerlines_main.cpp )
target_link_libraries( meshing_centerlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( centerlines.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_centerlines.cfg )

angiotk_add_component( meshing_remeshstl remeshstl_main.cpp )
target_link_libraries( meshing_remeshstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( remeshstlvmtk.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlvmtk.cfg )
configure_file( remeshstlgmsh.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_remeshstlgmsh.cfg )

angiotk_add_component( meshing_volumefromstlandcenterlines volumefromstlandcenterlines_main.cpp )
target_link_libraries( meshing_volumefromstlandcenterlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( volumefromstlandcenterlines.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstlandcenterlines.cfg )

angiotk_add_component( meshing_volumefromstl volumefromstl_main.cpp )
target_link_libraries( meshing_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
configure_file( volumefromstl.cfg ${AngioTk_BINARY_DIR}/cfg/meshing/test_volumefromstl.cfg )

angiotk_add_component( meshing_imagefromcenterlines imagefromcenterlines_main.cpp )
target_link_libraries( meshing_imagefromcenterlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_surfacefromimage surfacefromimage_main.cpp )
target_link_libraries( meshing_surfacefromimage ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_opensurface opensurface_main.cpp )
target_link_libraries( meshing_opensurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_smoothsurface smoothsurface_main.cpp )
target_link_libraries( meshing_smoothsurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_subdividesurface subdividesurface_main.cpp )
target_link_libraries( meshing_subdividesurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( meshing_centerlinesmanager centerlinesmanager_main.cpp )
target_link_libraries( meshing_centerlinesmanager ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)


ExternalData_expand_arguments( angioTkData testArgs DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation.stl} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation_remeshVMTKarea0.005.stl} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation_remeshVMTKarea0.005_centerlines.vtk} )
ExternalData_expand_arguments( angioTkData testArgs DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation.inletoutlet.desc} )
add_dependencies( meshing_volumefromstl angioTkData)
