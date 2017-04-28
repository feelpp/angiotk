
# Use external project module for building rorpo
include(ExternalProject)

# Get relative path from top source dir
string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_PATH_FROM_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})

# RORPO requires OpenMP, so we need to use gcc/g++
# looking for gcc/g++
execute_process(COMMAND which gcc OUTPUT_VARIABLE _GNU_C_COMPILER)
string(REPLACE "\n" "" _GNU_C_COMPILER ${_GNU_C_COMPILER})
execute_process(COMMAND which g++ OUTPUT_VARIABLE _GNU_CXX_COMPILER)
string(REPLACE "\n" "" _GNU_CXX_COMPILER ${_GNU_CXX_COMPILER})

# Checking version
# RORPO uses std::regex_iterator, which is not implemented before GCC 4.9
# Thus we have to ensure that we have the correct version
execute_process(COMMAND bash "-c" "`which gcc` --version" OUTPUT_VARIABLE _GNU_C_COMPILER_FULL_VERSION)
string(REGEX MATCH "[0-9]+[.][0-9]+[.][0-9]+" _GNU_C_COMPILER_DOTTED_VERSION ${_GNU_C_COMPILER_FULL_VERSION})
string(REPLACE "." "" _GNU_C_COMPILER_VERSION ${_GNU_C_COMPILER_DOTTED_VERSION})
execute_process(COMMAND bash "-c" "`which g++` --version" OUTPUT_VARIABLE _GNU_CXX_COMPILER_FULL_VERSION)
string(REGEX MATCH "[0-9]+[.][0-9]+[.][0-9]+" _GNU_CXX_COMPILER_DOTTED_VERSION ${_GNU_CXX_COMPILER_FULL_VERSION})
string(REPLACE "." "" _GNU_CXX_COMPILER_VERSION ${_GNU_CXX_COMPILER_DOTTED_VERSION})

if(_GNU_C_COMPILER_VERSION GREATER "490" AND _GNU_CXX_COMPILER_VERSION GREATER "490")
# Include RORPO as an external project
ExternalProject_Add(RORPO    # Name for custom target
  PREFIX "${CMAKE_BINARY_DIR}/Modules/Filtering/RORPO"
  GIT_REPOSITORY "https://github.com/path-openings/RORPO"
  CMAKE_ARGS -DCMAKE_C_COMPILER=${_GNU_C_COMPILER} -DCMAKE_CXX_COMPILER=${_GNU_CXX_COMPILER} -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}/${RELATIVE_PATH_FROM_TOP_SOURCE_DIR}
)
else()
    message(WARNING "RORPO uses std::regex_iterator which is not implemented properly before GCC 4.9. You are using gcc ${_GNU_C_COMPILER_VERSION} and g++ ${_GNU_CXX_COMPILER_VERSION}, thus it won't be built. Please update your alternatives for gcc and g++.")
endif()
