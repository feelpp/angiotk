#=========================================================
cmake_minimum_required (VERSION 3.1.3)

#=========================================================
#AngioTk = Virtual Angiography Toolkit
PROJECT(AngioTk)
#=========================================================
option(ANGIOTK_USE_SUPERBUILD "Build AngioTk and ITK/VMTK if required" ON )
if (ANGIOTK_USE_SUPERBUILD)
  set(ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION OFF)
else()
  set(ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION ON)
endif()
option(ANGIOTK_USE_SYSTEM_ITK "detect ITK on system" ${ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION})
option(ANGIOTK_USE_SYSTEM_VMTK "detect VMTK on system" ${ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION})
option(ANGIOTK_USE_SYSTEM_RORPO "detect RORPO on system" ${ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION})
if ( ANGIOTK_USE_SUPERBUILD )
  if ( NOT ANGIOTK_USE_SYSTEM_ITK )
    set(ANGIOTK_SUPERBUILD_ITK_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "ITK install prefix")
  endif()
  if ( NOT ANGIOTK_USE_SYSTEM_VMTK )
    set(ANGIOTK_SUPERBUILD_VMTK_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "VMTK install prefix")
  endif()
  if ( NOT ANGIOTK_USE_SYSTEM_RORPO )
    set(ANGIOTK_SUPERBUILD_RORPO_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "RORPO install prefix")
  endif()
endif()
OPTION(BUILD_EXAMPLES "Build examples" OFF)
OPTION(BUILD_DOXYGEN "Build Doxygen Documentation" OFF)
#=========================================================
# Shared libraries option
OPTION(BUILD_SHARED_LIBS "Build AngioTk with shared libraries." ON)
SET(AngioTk_BUILD_SHARED_LIBS ${BUILD_SHARED_LIBS})
#=========================================================
if (NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()
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
if ( 1 )
    if (ANGIOTK_USE_SUPERBUILD)
      set(FEELPP_DONT_SETUP_CMAKE 1)
    endif()
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

#=========================================================
# Find Python
if ( NOT PYTHONLIBS_FOUND )
  FIND_PACKAGE(PythonLibs)
endif()
if ( NOT PYTHONINTERP_FOUND )
  FIND_PACKAGE(PythonInterp)
endif()

#=========================================================
# Find ITK (required)
if ( NOT FEELPP_HAS_VTK )
  option(ANGIOTK_USE_SYSTEM_VTK "detect VTK on system" ${ANGIOTK_EXTERNAL_PACKAGES_USE_SYSTEM_DEFAULT_OPTION})
  if ( ANGIOTK_USE_SUPERBUILD AND (NOT ANGIOTK_USE_SYSTEM_VTK) )
    set(ANGIOTK_SUPERBUILD_VTK_INSTALL_PREFIX ${CMAKE_INSTALL_PREFIX} CACHE STRING "VTK install prefix")
  endif()
  if ( ANGIOTK_USE_SYSTEM_VTK )
    FIND_PACKAGE(VTK)
    IF(VTK_FOUND)
      message(STATUS "[VTK] Use file: ${VTK_USE_FILE}")
      INCLUDE(${VTK_USE_FILE})
    ELSE(VTK_FOUND)
      MESSAGE(FATAL_ERROR "Cannot build without VTK. Set VTK_DIR.")
    ENDIF(VTK_FOUND)
  endif()
endif()

#=========================================================
if ( ANGIOTK_USE_SYSTEM_VMTK )
# Find VMTK

# binary dir
unset(VMTK_EXECUTABLE_DIR)
FIND_PATH(VMTK_EXECUTABLE_DIR 
    NAMES vmtkcenterlines 
    #PATH
    HINTS
    ${VMTK_DIR}/bin
    $ENV{VMTK_DIR}/bin
    $ENV{VMTKHOME}/bin
    #NO_DEFAULT_PATH
    )

# include dir
FIND_PATH(VMTK_INCLUDE_DIR vtkvmtkConfigure.h 
  HINTS
  ${VMTK_DIR}/include/
  $ENV{VMTK_DIR}/include/
  $ENV{VMTKHOME}/include/
  PATH_SUFFIXES vmtk
  )

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

unset(VMTK_FOUND)
if ( VMTK_EXECUTABLE_DIR AND VMTK_INCLUDE_DIR AND VMTK_LIBRARIES )
  set( VMTK_FOUND 1)
  include_directories( ${VMTK_INCLUDE_DIR} )
  message(STATUS "found VMTK")
  message(STATUS "VMTK executables dir: ${VMTK_EXECUTABLE_DIR}")
  message(STATUS "VMTK included: ${VMTK_INCLUDE_DIR}")
  message(STATUS "VMTK libraries: ${VMTK_LIBRARIES}")
endif()

endif()

#=========================================================
# RORPO
if ( ANGIOTK_USE_SYSTEM_RORPO )

  unset(RORPO_EXECUTABLE_DIR)
  FIND_PATH(RORPO_EXECUTABLE_DIR
    NAMES RORPO_multiscale_usage
    HINTS
    ${RORPO_DIR}/bin
    )

  unset(RORPO_LIBRARIES)
  foreach( _rorpolib Nifti  RORPO  docopt)
    FIND_LIBRARY(RORPO_LIBRARY
      NAMES ${_rorpolib}
      HINTS
      ${RORPO_DIR}/lib
      )
    if ( RORPO_LIBRARY )
      if ( NOT RORPO_LIBRARIES )
	set(RORPO_LIBRARIES ${RORPO_LIBRARY})
      else()
	list( APPEND RORPO_LIBRARIES ${RORPO_LIBRARY} )
      endif()
      unset(RORPO_LIBRARY CACHE)
    endif()
  endforeach()

  unset(RORPO_FOUND)
  if ( RORPO_EXECUTABLE_DIR AND RORPO_LIBRARIES )
    set( RORPO_FOUND 1)
    message(STATUS "found RORPO")
    message(STATUS "RORPO executables dir: ${RORPO_EXECUTABLE_DIR}")
    message(STATUS "RORPO libraries: ${RORPO_LIBRARIES}")
  endif()
endif()

#=========================================================
if ( ANGIOTK_USE_SUPERBUILD )
  message(STATUS "AngioTk run SuperBuild")
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
  ADD_SUBDIRECTORY(Examples)
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
#add_subdirectory(Examples/MeshFromMRI)

include(CMake/AngioTkDirective.cmake)
