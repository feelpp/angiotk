
project(ResampleVolumesToBeIsotropic)

if ( ITK_FOUND )
  angiotk_add_component( angiotk_filtering_ResampleVolumesToBeIsotropic ResampleVolumesToBeIsotropic.cxx )
  target_link_libraries( angiotk_filtering_ResampleVolumesToBeIsotropic ${ITK_LIBRARIES})
endif()
