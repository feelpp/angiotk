
add_library( angiotk_meshing_volumefromstl
  SHARED
  volumefromstl.cpp
  AngioTkCenterlineField.cpp
  angiotkPolyBallModeller.cxx
  angiotkPolyBallLine.cxx
  centerlinesmanagerwindowinteractor.cpp
  centerlinesmanageriodata.cpp
  )
target_link_libraries( angiotk_meshing_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} ${VMTK_LIBRARIES} )
target_include_directories(angiotk_meshing_volumefromstl BEFORE PUBLIC ${VolumeFromSTL_SOURCE_DIR}/include/gmshHeadersMissing )
install(TARGETS angiotk_meshing_volumefromstl LIBRARY DESTINATION lib/ COMPONENT Libs)

if(APPLE)
  add_executable(angiotk_meshing_centerlinesmanagergui MACOSX_BUNDLE centerlinesmanagergui_main.cpp )
  target_link_libraries( angiotk_meshing_centerlinesmanagergui  ${VTK_LIBRARIES} angiotk_meshing_volumefromstl)
else()
  angiotk_add_component( angiotk_meshing_centerlinesmanagergui centerlinesmanagergui_main.cpp )
  #target_link_libraries( meshing_centerlinesmanagergui ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)
  target_link_libraries( angiotk_meshing_centerlinesmanagergui ${VTK_LIBRARIES} angiotk_meshing_volumefromstl)
endif()

angiotk_add_component( angiotk_meshing_createinletoutletdesc createinletoutletdesc_main.cpp )
target_link_libraries( angiotk_meshing_createinletoutletdesc ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_centerlines centerlines_main.cpp )
target_link_libraries( angiotk_meshing_centerlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_remeshstl remeshstl_main.cpp )
target_link_libraries( angiotk_meshing_remeshstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_volumefromstlandcenterlines volumefromstlandcenterlines_main.cpp )
target_link_libraries( angiotk_meshing_volumefromstlandcenterlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_pipeline_volumefromstl volumefromstl_main.cpp )
target_link_libraries( angiotk_meshing_pipeline_volumefromstl ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_imagefromcenterlines imagefromcenterlines_main.cpp )
target_link_libraries( angiotk_meshing_imagefromcenterlines ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_surfacefromimage surfacefromimage_main.cpp )
target_link_libraries( angiotk_meshing_surfacefromimage ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_imagesmanager imagesmanager_main.cpp )
target_link_libraries( angiotk_meshing_imagesmanager ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_opensurface opensurface_main.cpp )
target_link_libraries( angiotk_meshing_opensurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_smoothsurface smoothsurface_main.cpp )
target_link_libraries( angiotk_meshing_smoothsurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_subdividesurface subdividesurface_main.cpp )
target_link_libraries( angiotk_meshing_subdividesurface ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_centerlinesmanager centerlinesmanager_main.cpp )
target_link_libraries( angiotk_meshing_centerlinesmanager ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)

angiotk_add_component( angiotk_meshing_tubularextension tubularextension_main.cpp )
target_link_libraries( angiotk_meshing_tubularextension ${FEELPP_LIBRARY} ${FEELPP_LIBRARIES} angiotk_meshing_volumefromstl)


