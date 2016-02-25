
include(ExternalProject)

# Force C++ standard to be 98
# This sample will not build with c++11 or c++14
# Useful if you want to directly built the sample inside source
#set(_PREV_CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++98")

# If ITK is available, we fetch the Resamplevolumestobeisotropic example from ITK
# This is done to overcome the fact that the examples might not be built in 
# the installed ITK version
if("${ITK_FOUND}")
    ExternalProject_Add(ResampleVolumesToBeIsotropic
        PREFIX "${CMAKE_BINARY_DIR}"
        DOWNLOAD_COMMAND 
            mkdir -p ${CMAKE_BINARY_DIR}/src/ResampleVolumesToBeIsotropic/ && cd ${CMAKE_BINARY_DIR}/src/ResampleVolumesToBeIsotropic/
            && wget https://raw.githubusercontent.com/Kitware/ITK/master/Examples/Filtering/ResampleVolumesToBeIsotropic.cxx -O ResampleVolumesToBeIsotropic.cxx
            && echo "cmake_minimum_required(VERSION 2.8)" > CMakeLists.txt
            && echo "find_package(ITK)" >> CMakeLists.txt
            && echo "include(${ITK_USE_FILE})" >> CMakeLists.txt
            && echo "add_executable(ResampleVolumesToBeIsotropic ResampleVolumesToBeIsotropic.cxx )" >> CMakeLists.txt
            && echo "target_link_libraries(ResampleVolumesToBeIsotropic ${ITK_LIBRARIES})" >> CMakeLists.txt
            && echo "install(TARGETS ResampleVolumesToBeIsotropic RUNTIME DESTINATION Modules/Filtering/ITK/bin)" >> CMakeLists.txt
        CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${CMAKE_INSTALL_PREFIX}
        )
else()
    message("ITK was not found on your system. Isotropic resampling is thus disabled.")
endif()

# Restore previous flags
#set(${CMAKE_CXX_FLAGS} _PREV_CMAKE_CXX_FLAGS)
