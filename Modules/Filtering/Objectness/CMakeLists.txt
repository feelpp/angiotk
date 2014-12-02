project(Filtering)

include_directories(${Filtering_SOURCE_DIR}/include)

add_subdirectory( src )

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
