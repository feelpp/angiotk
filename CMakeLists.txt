#=========================================================
CMAKE_MINIMUM_REQUIRED(VERSION 2.8.7 FATAL_ERROR)

#=========================================================
#AngioTk = Virtual Angiography Toolkit
PROJECT(AngioTk)
#=========================================================
option(ANGIOTK_USE_SYSTEM_ITK "detect ITK on system" OFF)
option(ANGIOTK_USE_SYSTEM_VMTK "detect VMTK on system" OFF)
option(ANGIOTK_USE_SUPERBUILD "Build AngioTk and ITK/VMTK if required" ON )
OPTION(BUILD_EXAMPLES "Build examples" OFF)
OPTION(BUILD_DOXYGEN "Build Doxygen Documentation" OFF)
#=========================================================
# Shared libraries option
OPTION(BUILD_SHARED_LIBS "Build AngioTk with shared libraries." ON)
SET(AngioTk_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
#=========================================================

#=========================================================
#Support for the CTest dashboard testing system
INCLUDE(CTest)
#=========================================================

IF(BUILD_TESTING)
  ENABLE_TESTING()
  INCLUDE(CMake/AngioTkExternalData.cmake)
  
  # Sets the data root
  set(AngiotTk_DATA_ROOT ${AngioTk_SOURCE_DIR}/Testing/Data)
  mark_as_advanced(AngiotTk_DATA_ROOT)
ENDIF(BUILD_TESTING)
#=========================================================

#=========================================================
SET(AngioTk_VERSION_MAJOR "0")
SET(AngioTk_VERSION_MINOR "1")
SET(AngioTk_VERSION_PATCH "0")

#=========================================================

#=========================================================
# Output directories.
IF(NOT LIBRARY_OUTPUT_PATH)
  SET (LIBRARY_OUTPUT_PATH ${AngioTk_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all libraries.")
ENDIF(NOT LIBRARY_OUTPUT_PATH)
IF(NOT EXECUTABLE_OUTPUT_PATH)
  SET (EXECUTABLE_OUTPUT_PATH ${AngioTk_BINARY_DIR}/bin CACHE INTERNAL "Single output directory for building all executables.")
ENDIF(NOT EXECUTABLE_OUTPUT_PATH)
MARK_AS_ADVANCED(LIBRARY_OUTPUT_PATH EXECUTABLE_OUTPUT_PATH)

SET(AngioTk_LIBRARY_PATH "${LIBRARY_OUTPUT_PATH}")
SET(AngioTk_EXECUTABLE_PATH "${EXECUTABLE_OUTPUT_PATH}")
set(AngioTk_TEST_OUTPUT_DIR "${AngioTk_BINARY_DIR}/Testing/Temporary")
SET(CXX_TEST_PATH ${EXECUTABLE_OUTPUT_PATH})

#=========================================================
# Remove some MS Visual c++ flags
IF(MSVC)
  ADD_DEFINITIONS(-D_CRT_SECURE_NO_WARNINGS -D_CRT_NONSTDC_NO_DEPRECATE -D_SCL_SECURE_NO_WARNINGS)
ENDIF(MSVC)

# Flag required for shared libraries
IF ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
  SET (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fPIC")
  SET (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC")
ENDIF ("${CMAKE_SYSTEM_PROCESSOR}" STREQUAL "x86_64")
#=========================================================

#=========================================================
# Find ITK (required)
if ( ANGIOTK_USE_SYSTEM_ITK )
FIND_PACKAGE(ITK)
IF(ITK_FOUND)
  message(STATUS "[ITK] Use file: ${ITK_USE_FILE}")
  INCLUDE(${ITK_USE_FILE})
  if(${ITK_VERSION_MAJOR} EQUAL 3)
    message(FATAL_ERROR "AngioTk only works with ITK version 4 and higher")
  endif(${ITK_VERSION_MAJOR} EQUAL 3)
ELSE(ITK_FOUND)
  MESSAGE(FATAL_ERROR "Cannot build without ITK. Set ITK_DIR.")
ENDIF(ITK_FOUND)
endif()

#=========================================================
# Find Feel++
OPTION(ANGIOTK_USE_OLD_FEELPP_CMAKE_DETECTION "use old cmake feel++ detection" OFF)
if ( ANGIOTK_USE_OLD_FEELPP_CMAKE_DETECTION )
  FIND_PATH(FEELPP_CMAKE_MODULES FindFeel++.cmake PATH $ENV{FEELPP_DIR}/share/feel/cmake/modules/   /usr/share/feel/cmake/modules/ /usr/local/share/feel/cmake/modules/ )
  if ( FEELPP_CMAKE_MODULES )
    set(CMAKE_MODULE_PATH ${FEELPP_CMAKE_MODULES})
    Find_Package(Feel++)
    set(Feel++_FOUND 1)
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g0 -O3 -DNDEBUG  -finline-functions -pedantic  -ftemplate-depth-256 -Wno-inline" )

    if(BUILD_MODULE_Meshing AND (NOT ANN_INCLUDE_DIR OR NOT ANN_LIB))
      message(FATAL_ERROR "ANN does not seem to be installed on this platform. Please install Feel++ with ANN support.")
    endif()
  else()
    # The meshing module uses Feel++
    # If it is not installed, make cmake end up with an error
    if(BUILD_MODULE_Meshing)
      message(FATAL_ERROR "Feel++ does not seem to be installed on this platform. Please set the FEELPP_DIR variable.")
    endif()
  endif()
else()
  if ( ${CMAKE_SOURCE_DIR} STREQUAL ${CMAKE_CURRENT_SOURCE_DIR} )
    find_package(Feel++ PATHS
      ${FEELPP_DIR}/share/feelpp/feel/cmake/modules
      $ENV{FEELPP_DIR}/share/feelpp/feel/cmake/modules
      /usr/share/feelpp/feel/cmake/modules /usr/local/share/feelpp/feel/cmake/modules /opt/share/feelpp/feel/cmake/modules)
    if(NOT FEELPP_FOUND)
      message(FATAL_ERROR "Feel++ was not found on your system. Make sure to install it and specify the FEELPP_DIR to reference the installation directory.") 
    endif()
    if(BUILD_MODULE_Meshing AND (NOT FEELPP_HAS_ANN))
      message(FATAL_ERROR "ANN does not seem to be installed on this platform. Please install Feel++ with ANN support.")
    endif()

    if ( 0 )
    # Check that the COMPILER flags are compatible, in case we don't compile with the same compiler as Feel++ was built with
    # Clang and GCC might have flag differences
    include(CheckCXXCompilerFlag)
    # store the previous value of CXX flags
    string(REGEX REPLACE "[ ]+" " " _CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
    set(CMAKE_CXX_FLAGS "")
    string(REPLACE " " ";" _CXX_FLAGS ${_CMAKE_CXX_FLAGS})
    set(REMOVED_FLAGS "")

    # Remove flags if needed
    set(_counter 0)
    foreach(flag ${_CXX_FLAGS})
        CHECK_CXX_COMPILER_FLAG("${flag}" FLAG_${_counter})
        if(${FLAG_${_counter}})
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}")
        else()
            set(REMOVED_FLAGS "${REMOVED_FLAGS} ${flag}")
        endif()
        math(EXPR _counter "${_counter} + 1")
    endforeach()

    # Inform the user about potentially removed flags
    if(REMOVED_FLAGS)
        message(WARNING "Initial C++ flags: ${_CMAKE_CXX_FLAGS}")
        message(WARNING "Removed flags: ${REMOVED_FLAGS}")
    endif()
    unset(_CMAKE_CXX_FLAGS)
    unset(REMOVED_FLAGS)
    endif()
endif()
endif()

#=========================================================
# Find Python
#FIND_PACKAGE(PythonLibs)
#FIND_PACKAGE(PythonInterp)
# Python 3 not supported by vmtk <=1.3
FIND_PACKAGE(PythonLibs 2)
FIND_PACKAGE(PythonInterp 2)

#=========================================================
if ( ANGIOTK_USE_SYSTEM_VMTK )
# Find VMTK

# binary dir
FIND_PATH(VMTK_BINARY_DIR 
    NAMES vmtkcenterlines 
    PATH 
    ${VMTK_DIR}/bin
    $ENV{VMTK_DIR}/bin
    $ENV{VMTKHOME}/bin
    )
# if ( NOT VMTK_BINARY_DIR )
#    message(FATAL_ERROR "VMTK_BINARY_DIR not find")
# endif()

# include dir
FIND_PATH(VMTK_INCLUDE_DIR vtkvmtkConfigure.h 
  HINTS
  ${VMTK_DIR}/include/
  $ENV{VMTK_DIR}/include/
  $ENV{VMTKHOME}/include/
  PATH_SUFFIXES vmtk
  )
#FIND_PATH(VMTK_INCLUDE_DIR vmtk/vtkvmtkConfigure.h HINTS $ENV{VMTK_DIR}/include/  )
# if ( NOT VMTK_INCLUDE_DIR )
#   message(FATAL_ERROR "VMTK_INCLUDE_DIR not find")
# endif()

# libraries path
unset(VMTK_LIBRARIES)
foreach(_vmtklib vtkvmtkCommon vtkvmtkComputationalGeometry)
  unset(VMTK_LIBRARY CACHE)
  FIND_LIBRARY(VMTK_LIBRARY
      NAMES ${_vmtklib} 
      HINTS
      ${VMTK_DIR}/lib/
      $ENV{VMTK_DIR}/lib/
      $ENV{VMTKHOME}/lib/
      PATH_SUFFIXES vmtk
      )
  if ( VMTK_LIBRARY )
    set( VMTK_LIBRARIES ${VMTK_LIBRARIES} ${VMTK_LIBRARY} )
  endif()
endforeach()

if ( VMTK_BINARY_DIR AND VMTK_INCLUDE_DIR AND VMTK_LIBRARIES )
  set( VMTK_FOUND 1)
else()
  set( VMTK_FOUND 0)
endif()

if ( ${VMTK_FOUND} )
  include_directories( ${VMTK_INCLUDE_DIR} )
  message(STATUS "VMTK binary dir: ${VMTK_BINARY_DIR}")
  message(STATUS "VMTK included: ${VMTK_INCLUDE_DIR}")
  message(STATUS "VMTK libraries: ${VMTK_LIBRARIES}")
endif()

endif( ANGIOTK_USE_SYSTEM_VMTK )

#=========================================================
if ( ANGIOTK_USE_SUPERBUILD )
  include(${AngioTk_SOURCE_DIR}/CMake/AngioTkSuperBuild.cmake)
  return()
endif()


#=========================================================
# some checks
if(BUILD_MODULE_Meshing)
  if ( NOT VMTK_FOUND )
    message(FATAL_ERROR "VMTK does not seem to be installed on this platform.
    Please set the VMTK_DIR variable or use the vmtk_env.sh script provided
    with VMTK and export VMTKHOME.")
  endif()
endif()



#=========================================================
# AngioTk uses KWStyle for checking the coding style
INCLUDE(${AngioTk_SOURCE_DIR}/CMake/KWStyle/KWStyle.cmake)
#=========================================================


#=========================================================
# If choose to build documentation, then search for Doxygen executables.

IF(BUILD_DOXYGEN)
  ADD_SUBDIRECTORY(documentation/Doxygen)
ENDIF(BUILD_DOXYGEN)
#=========================================================

#=========================================================
# If choose to build examples, then look into Examples folder.
IF(BUILD_EXAMPLES)
  ADD_SUBDIRECTORY(Examples/Python)
ENDIF(BUILD_EXAMPLES)
#=========================================================

#=========================================================
# Installation variables
if(NOT AngioTK_INSTALL_RUNTIME_DIR)
  set(AngioTK_INSTALL_RUNTIME_DIR bin)
endif()
if(NOT AngioTK_INSTALL_LIB_DIR)
  set(AngioTK_INSTALL_LIB_DIR lib)
endif()
if(NOT AngioTK_INSTALL_ARCHIVE_DIR)
  set(AngioTK_INSTALL_ARCHIVE_DIR lib)
endif()
if(NOT AngioTK_INSTALL_INCLUDE_DIR)
  set(AngioTK_INSTALL_INCLUDE_DIR include/AngioTk-${AngioTK_VERSION_MAJOR}.${AngioTK_VERSION_MINOR})
endif()
if(NOT AngioTK_INSTALL_PACKAGE_DIR)
  set(AngioTK_INSTALL_PACKAGE_DIR "lib/cmake/AngioTk-${AngioTK_VERSION_MAJOR}.${AngioTK_VERSION_MINOR}")
endif()
#=========================================================

INCLUDE_DIRECTORIES(${AngioTK_INCLUDE_DIRS})

# Useful generic macros for modules configuration  
INCLUDE("${AngioTk_SOURCE_DIR}/CMake/AngioTkMacros.cmake")

# Initiate modules configuration
INCLUDE("${AngioTk_SOURCE_DIR}/CMake/AngioTkModules.cmake")

ExternalData_Add_Target(angioTkData)

add_subdirectory(Data)
add_subdirectory(Pipeline)
add_subdirectory(Examples/MeshFromMRI)

include(CMake/AngioTkDirective.cmake)
