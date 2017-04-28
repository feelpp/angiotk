
# Use external project module for building rorpo
include(ExternalProject)

# Get relative path from top source dir
string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_PATH_FROM_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# RORPO requires OpenMP, so we need to use gcc/g++
execute_process(COMMAND which gcc OUTPUT_VARIABLE _GNU_C_COMPILER)
string(REPLACE "\n" "" _GNU_C_COMPILER ${_GNU_C_COMPILER})
execute_process(COMMAND which g++ OUTPUT_VARIABLE _GNU_CXX_COMPILER)
string(REPLACE "\n" "" _GNU_CXX_COMPILER ${_GNU_CXX_COMPILER})

# Include RORPO as an external project
ExternalProject_Add(RORPO    # Name for custom target
  PREFIX "${CMAKE_BINARY_DIR}/Modules/Filtering/RORPO"
  GIT_REPOSITORY "https://github.com/path-openings/RORPO"
  CMAKE_ARGS -DCMAKE_C_COMPILER=${_GNU_C_COMPILER} -DCMAKE_CXX_COMPILER=${_GNU_CXX_COMPILER} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/${RELATIVE_PATH_FROM_TOP_SOURCE_DIR}
)
