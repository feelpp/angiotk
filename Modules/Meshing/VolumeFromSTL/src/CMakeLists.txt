
set(ANGIOTK_MODULE_MESHING_EXTERNAL_LIBRARIES)
if ( FEELPP_FOUND )
  list(APPEND ANGIOTK_MODULE_MESHING_EXTERNAL_LIBRARIES ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES})
endif()
if ( VMTK_FOUND )
  list(APPEND ANGIOTK_MODULE_MESHING_EXTERNAL_LIBRARIES ${VMTK_LIBRARIES})
endif()
if ( VTK_FOUND )
  list(APPEND ANGIOTK_MODULE_MESHING_EXTERNAL_LIBRARIES ${VTK_LIBRARIES})
endif()

add_library( angiotk_meshing_volumefromstl
  SHARED
  volumefromstl.cpp
  AngioTkCenterlineField.cpp
  angiotkPolyBallModeller.cxx
  angiotkPolyBallLine.cxx
  centerlinesmanagerwindowinteractor.cpp
  centerlinesmanageriodata.cpp
  )
target_link_libraries( angiotk_meshing_volumefromstl ${ANGIOTK_MODULE_MESHING_EXTERNAL_LIBRARIES}) # ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${VMTK_LIBRARIES} )
target_include_directories(angiotk_meshing_volumefromstl BEFORE PUBLIC ${VolumeFromSTL_SOURCE_DIR}/include/gmshHeadersMissing )
install(TARGETS angiotk_meshing_volumefromstl LIBRARY DESTINATION lib/ COMPONENT Libs)

set(ANGIOTK_MODULE_MESHING_LIBRARIES angiotk_meshing_volumefromstl ) # ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

if(APPLE)
  add_executable(angiotk_meshing_centerlinesmanagergui MACOSX_BUNDLE centerlinesmanagergui_main.cpp )
  target_link_libraries( angiotk_meshing_centerlinesmanagergui ${ANGIOTK_MODULE_MESHING_LIBRARIES} )
else()
  angiotk_add_component( angiotk_meshing_centerlinesmanagergui centerlinesmanagergui_main.cpp )
  target_link_libraries( angiotk_meshing_centerlinesmanagergui ${ANGIOTK_MODULE_MESHING_LIBRARIES})
endif()

angiotk_add_component( angiotk_meshing_createinletoutletdesc createinletoutletdesc_main.cpp )
target_link_libraries( angiotk_meshing_createinletoutletdesc ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_centerlines centerlines_main.cpp )
target_link_libraries( angiotk_meshing_centerlines ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_remeshstl remeshstl_main.cpp )
target_link_libraries( angiotk_meshing_remeshstl ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_volumefromstlandcenterlines volumefromstlandcenterlines_main.cpp )
target_link_libraries( angiotk_meshing_volumefromstlandcenterlines ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_pipeline_volumefromstl volumefromstl_main.cpp )
target_link_libraries( angiotk_meshing_pipeline_volumefromstl ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_imagefromcenterlines imagefromcenterlines_main.cpp )
target_link_libraries( angiotk_meshing_imagefromcenterlines ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_surfacefromimage surfacefromimage_main.cpp )
target_link_libraries( angiotk_meshing_surfacefromimage ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_imagesmanager imagesmanager_main.cpp )
target_link_libraries( angiotk_meshing_imagesmanager ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_opensurface opensurface_main.cpp )
target_link_libraries( angiotk_meshing_opensurface ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_smoothsurface smoothsurface_main.cpp )
target_link_libraries( angiotk_meshing_smoothsurface ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_subdividesurface subdividesurface_main.cpp )
target_link_libraries( angiotk_meshing_subdividesurface ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_centerlinesmanager centerlinesmanager_main.cpp )
target_link_libraries( angiotk_meshing_centerlinesmanager ${ANGIOTK_MODULE_MESHING_LIBRARIES} )

angiotk_add_component( angiotk_meshing_tubularextension tubularextension_main.cpp )
target_link_libraries( angiotk_meshing_tubularextension ${ANGIOTK_MODULE_MESHING_LIBRARIES} )


