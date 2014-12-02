# Module Meshing
# Label volume to Mesh for CFD

project(Meshing)
SET( BUILD_MODULE_Meshing OFF )

include_directories(${Meshing_SOURCE_DIR}/include)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()

