# Module Meshing
# Label volume to Mesh for CFD

project(RemeshBasedOnCenterlines)

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

FIND_PACKAGE(PythonLibs)
FIND_PACKAGE(PythonInterp)
configure_file(angiotkMeshingConfig.h.in angiotkMeshingConfig.h)
include_directories(${RemeshBasedOnCenterlines_BINARY_DIR})

include_directories(${RemeshBasedOnCenterlines_SOURCE_DIR}/include)



if(BUILD_TESTING)
  add_subdirectory(test)
endif()

endif( Feel++_FOUND )


