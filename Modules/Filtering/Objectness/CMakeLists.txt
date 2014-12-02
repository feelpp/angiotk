project(Filtering)

include_directories(${Filtering_SOURCE_DIR}/include)

if(BUILD_TESTING)
  add_subdirectory(test)
endif()
