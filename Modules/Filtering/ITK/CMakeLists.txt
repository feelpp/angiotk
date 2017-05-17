
project(ResampleVolumesToBeIsotropic)

if ( ITK_FOUND )
  angiotk_add_component( filtering_ResampleVolumesToBeIsotropic ResampleVolumesToBeIsotropic.cxx )
  target_link_libraries( filtering_ResampleVolumesToBeIsotropic ${ITK_LIBRARIES})
endif()
