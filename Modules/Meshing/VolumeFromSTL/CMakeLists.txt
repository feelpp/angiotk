# Module Meshing
# Label volume to Mesh for CFD

project(VolumeFromSTL)

#set(ENV{FEELPP_DIR} /Users/vincentchabannes/github-experimental/feelpp-install/)
FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH $ENV{FEELPP_DIR}/share/feel/cmake/modules/   /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
if ( FEELPP_CMAKE_MODULES )
  set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
  Find_Package(Feel++)
  set(Feel++_FOUND 1)
else()
  #message(FATAL_ERROR "Feel++ does not seem to have been installed on this platform")
endif()


if ( Feel++_FOUND )

##################################################
# Python
FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(PythonInterp)

##################################################
# VMTK
#set(ENV{VMTK_DIR} /Users/vincentchabannes/packages/vmtk/vmtk.build2/Install/bin/)
FIND_PATH(VMTK_BINARY_DIR NAMES vmtkcenterlines PATH $ENV{VMTK_DIR}/bin )
if ( NOT VMTK_BINARY_DIR )
  message(FATAL_ERROR "VMTK_BINARY_DIR not find")
endif()

configure_file(angiotkMeshingConfig.h.in angiotkMeshingConfig.h)
include_directories(${VolumeFromSTL_BINARY_DIR})

include_directories(${VolumeFromSTL_SOURCE_DIR}/include)


add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

endif( Feel++_FOUND )


