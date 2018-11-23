
include(ExternalProject)


set(ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS)
set(ANGIOTK_EXTERNALPROJECT_DEPENDS)

set(ANGIOTK_SUPERBUILD_SETUP_TEXT)

if (APPLE)
  set( ANGIOTK_LD_LIBRARY_PATH_VAR_NAME "DYLD_LIBRARY_PATH")
else()
  set( ANGIOTK_LD_LIBRARY_PATH_VAR_NAME "LD_LIBRARY_PATH")
endif()
#####################
# ITK
#####################
if ( NOT ANGIOTK_USE_SYSTEM_ITK )
  ExternalProject_Add(AngioTk_ExternalPackages_ITK    # Name for custom target
    PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/ITK"
    GIT_REPOSITORY "https://github.com/InsightSoftwareConsortium/ITK.git"
    GIT_TAG v4.11.1
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
    -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DModule_ITKReview=ON -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/ITK
    #UPDATE_DISCONNECTED 1
    UPDATE_COMMAND ""
    )

  set( ITK_INSTALL_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/ITK)
  set( ITK_DIR ${ITK_INSTALL_DIR}/lib/cmake/ITK-4.11)
  set( ITK_FOUND 1)
  set( ANGIOTK_HAS_ITK_FROM_SUPERBUILD 1 )
  ExternalProject_Get_Property(AngioTk_ExternalPackages_ITK  BINARY_DIR )
  set( ANGIOTK_SUPERBUILD_ITK_BINARY_DIR ${BINARY_DIR} )
  list(APPEND ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_ITK)
  list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_ITK)

  set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export ${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}=${ITK_INSTALL_DIR}/lib:$${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}\n")
endif()

#####################
# VTK
#####################
if ( (NOT ANGIOTK_USE_SYSTEM_VTK) AND (NOT FEELPP_HAS_VTK) )
  set(VTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
     #-DVTK_USE_PARALLEL=ON -DVTK_USE_MPI=ON
     -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF
     -DVTK_WRAP_PYTHON=ON -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE} -DPYTHON_INCLUDE_DIR=${PYTHON_INCLUDE_DIR} -DPYTHON_LIBRARY=${PYTHON_LIBRARY}
     -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/VTK
     )
  ExternalProject_Add(AngioTk_ExternalPackages_VTK    # Name for custom target
    PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/VTK"
    GIT_REPOSITORY "https://github.com/Kitware/VTK.git"
    GIT_TAG v7.1.0
    CMAKE_ARGS ${VTK_SUPERBUILD_CMAKE_ARGS}
    #UPDATE_DISCONNECTED 1
    UPDATE_COMMAND ""
    )

  set( VTK_INSTALL_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/VTK)
  set( VTK_DIR ${VTK_INSTALL_DIR}/lib/cmake/vtk-7.1)
  set( VTK_FOUND 1)
  set( ANGIOTK_HAS_VTK_FROM_SUPERBUILD 1 )
  ExternalProject_Get_Property(AngioTk_ExternalPackages_VTK  BINARY_DIR )
  set( ANGIOTK_SUPERBUILD_VTK_BINARY_DIR ${BINARY_DIR} )
  list(APPEND ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_VTK)
  list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_VTK)

  set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export ${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}=${VTK_INSTALL_DIR}/lib:$${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}\n")
  set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export PYTHONPATH=${VTK_INSTALL_DIR}/lib/python2.7/site-packages:$PYTHONPATH\n")
endif()
#####################
# VMTK
#####################
if ( NOT ANGIOTK_USE_SYSTEM_VMTK )
  set( VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
  list( APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DVMTK_USE_SUPERBUILD=OFF)
  list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_ITK=ON -DITK_DIR=${ITK_DIR} )
  if ( FEELPP_HAS_PARAVIEW )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_VTK=ON -DParaView_DIR=${FEELPP_PARAVIEW_DIR} -DVTK_DIR=${FEELPP_PARAVIEW_DIR} )
  elseif( FEELPP_HAS_VTK )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_VTK=ON -DVTK_DIR=${FEELPP_VTK_DIR} )
  elseif( VTK_FOUND )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_VTK=ON -DVTK_DIR=${VTK_DIR} )
  endif()
  if ( PYTHON_EXECUTABLE )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DPYTHON_EXECUTABLE=${PYTHON_EXECUTABLE})
  endif()
  IF (APPLE)
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_MACOSX_RPATH=1)
  endif()
  list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/VMTK)

  ExternalProject_Add(AngioTk_ExternalPackages_VMTK    # Name for custom target
    PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/VMTK"
    GIT_REPOSITORY "https://github.com/vmtk/vmtk.git"
    GIT_TAG 43fab892721ed540df36f088f25765c3ef4df3c2 #v1.4.0 # v1.3.2
    PATCH_COMMAND git apply ${AngioTk_SOURCE_DIR}/CMake/AngioTk_VMTK_v1.3.2.patch ${AngioTk_SOURCE_DIR}/CMake/AngioTk_VMTK_vtkVersionLess8.patch
    CMAKE_ARGS ${VMTK_SUPERBUILD_CMAKE_ARGS}
    #UPDATE_DISCONNECTED 1
    UPDATE_COMMAND ""
    DEPENDS ${ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS}
    )

  set(VMTK_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/VMTK)
  set( VMTK_FOUND 1)
  set( ANGIOTK_HAS_VMTK_FROM_SUPERBUILD 1 )
  ExternalProject_Get_Property(AngioTk_ExternalPackages_VMTK  BINARY_DIR )
  set( ANGIOTK_SUPERBUILD_VMTK_BINARY_DIR ${BINARY_DIR} )
  list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_VMTK)

  set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export ${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}=${VMTK_DIR}/lib:$${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}\n")
  set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export PYTHONPATH=${VMTK_DIR}/lib/python2.7/site-packages:$PYTHONPATH\n")
endif()

######################
# RORPO
#####################

if ( NOT ANGIOTK_USE_SYSTEM_RORPO )

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

    set( RORPO_SUPERBUILD_CMAKE_ARGS -DCMAKE_CXX_COMPILER=${_GNU_CXX_COMPILER} -DCMAKE_C_COMPILER=${_GNU_C_COMPILER}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE})
    list(APPEND RORPO_SUPERBUILD_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/RORPO)
    ExternalProject_Add(AngioTk_ExternalPackages_RORPO
      PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/RORPO"
      GIT_REPOSITORY "https://github.com/path-openings/RORPO"
      GIT_TAG "v1.0"
      CMAKE_ARGS ${RORPO_SUPERBUILD_CMAKE_ARGS}
      UPDATE_COMMAND ""
      )

    set( RORPO_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/RORPO )
    set( RORPO_FOUND 1)
    set( ANGIOTK_HAS_RORPO_FROM_SUPERBUILD 1 )
    ExternalProject_Get_Property(AngioTk_ExternalPackages_RORPO  BINARY_DIR )
    set( ANGIOTK_SUPERBUILD_RORPO_BINARY_DIR ${BINARY_DIR} )
    list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_RORPO)

    set(ANGIOTK_SUPERBUILD_SETUP_TEXT "${ANGIOTK_SUPERBUILD_SETUP_TEXT} export ${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}=${RORPO_DIR}/lib:$${ANGIOTK_LD_LIBRARY_PATH_VAR_NAME}\n")
  endif()
endif()


#####################
# AngioTk
#####################

set( ANGIOTK_SUPERBUILD_CMAKE_ARGS
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER} -DCMAKE_CXX_FLAGS=${CMAKE_CXX_FLAGS}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}
  -DBUILD_EXAMPLES=${BUILD_EXAMPLES} -DBUILD_DOXYGEN=${BUILD_DOXYGEN}
  -DANGIOTK_USE_SUPERBUILD=OFF
  -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/AngioTkSuperBuild/install
  )
if ( FEELPP_FOUND )
  list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DFEELPP_DIR=${FEELPP_DIR})
endif()
if ( ITK_FOUND )
  list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DANGIOTK_USE_SYSTEM_ITK=ON -DITK_DIR=${ITK_DIR})
endif()
if ( VMTK_FOUND )
  list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DANGIOTK_USE_SYSTEM_VMTK=ON -DVMTK_DIR=${VMTK_DIR})
endif()
if ( RORPO_FOUND )
  list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DANGIOTK_USE_SYSTEM_RORPO=ON  -DRORPO_DIR=${RORPO_DIR})
endif()
if ( VTK_FOUND )
  list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DANGIOTK_USE_SYSTEM_VTK=ON -DVTK_DIR=${VTK_DIR})
endif()

file( GLOB modules ${AngioTk_SOURCE_DIR}/Modules/* ) 
foreach( module ${modules} )
  if( IS_DIRECTORY ${module} )
    get_filename_component( moduleName ${module} NAME)
    if( ${BUILD_MODULE_${moduleName}} )
      list(APPEND ANGIOTK_SUPERBUILD_CMAKE_ARGS -DBUILD_MODULE_${moduleName}=ON)
    endif()
  endif()
endforeach()

ExternalProject_Add(AngioTkSuperBuild    # Name for custom target
  PREFIX "${AngioTk_BINARY_DIR}/AngioTkSuperBuild/build"
  SOURCE_DIR "${AngioTk_SOURCE_DIR}"
  DOWNLOAD_COMMAND ""
  CMAKE_ARGS ${ANGIOTK_SUPERBUILD_CMAKE_ARGS}
  #UPDATE_DISCONNECTED 1
  #UPDATE_COMMAND ""
  DEPENDS ${ANGIOTK_EXTERNALPROJECT_DEPENDS}
  #STEP_TARGETS configure build
  )

set( ANGIOTK_HAS_ANGIOTK_FROM_SUPERBUILD 1 )
ExternalProject_Get_Property(AngioTkSuperBuild  BINARY_DIR )
set( ANGIOTK_SUPERBUILD_ANGIOTK_BINARY_DIR ${BINARY_DIR} )

configure_file(${AngioTk_SOURCE_DIR}/CMake/AngioTkSuperBuildInstall.cmake.in ${AngioTk_BINARY_DIR}/CMake/AngioTkSuperBuildInstall.cmake  @ONLY)
install(SCRIPT ${AngioTk_BINARY_DIR}/CMake/AngioTkSuperBuildInstall.cmake)

file(WRITE  ${ANGIOTK_SUPERBUILD_ANGIOTK_BINARY_DIR}/AngioTkSetupEnv.sh ${ANGIOTK_SUPERBUILD_SETUP_TEXT})
