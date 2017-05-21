# Module Meshing
# Label volume to Mesh for CFD

project(VolumeFromSTL)

if ( Feel++_FOUND AND VMTK_FOUND AND PYTHONINTERP_FOUND )

  include_directories(${VolumeFromSTL_BINARY_DIR})

  include_directories(${VolumeFromSTL_SOURCE_DIR}/include)

  add_subdirectory(src)

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()

endif()


