project(component_tree_segmentation)
#cmake_minimum_required(VERSION 2.8)
set(CMAKE_CXX_FLAGS "-std=c++11")

message( STATUS "Reading CT_Segmentation CMakeFiles.txt" )

#include_directories(${CT_Segmentation_BINARY_DIR})
#include_directories(${CT_Segmentation_SOURCE_DIR})

# Add files
#file( GLOB_RECURSE HEADER_FILES *.h *.hpp)
#file( GLOB_RECURSE SOURCE_FILES *.c *.cpp )

#message( STATUS "Headers: ${HEADER_FILES}" )
#message( STATUS "Sources: ${SOURCE_FILES}" )

#add_executable(${PROJECT_NAME} ${SOURCE_FILES} ${HEADER_FILES})

include_directories( nifti )
include_directories( include )

add_library( cts_nifti
  SHARED
  nifti/nifti1_io.c
  nifti/znzlib.c
  )
install(TARGETS cts_nifti
    LIBRARY DESTINATION Modules/CT_Segmentation/lib
)
#target_link_libraries( cts_nifti )


angiotk_add_component( ${PROJECT_NAME} main.cpp )
target_link_libraries( ${PROJECT_NAME} cts_nifti )

angiotk_add_component( stloff stloff_main.cpp )