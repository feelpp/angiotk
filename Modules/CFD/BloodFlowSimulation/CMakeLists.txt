# Module Meshing
# Label volume to Mesh for CFD

project(BloodFlowSimulation)

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

SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O3 -DNDEBUG  -finline-functions -pedantic  -ftemplate-depth-256 -Wno-inline" )

include_directories(${BloodFlowSimulation_SOURCE_DIR}/include)

add_subdirectory(src)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

endif( Feel++_FOUND )


