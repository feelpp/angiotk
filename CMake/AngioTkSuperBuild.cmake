
include(ExternalProject)


set(ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS)
set(ANGIOTK_EXTERNALPROJECT_DEPENDS)

#####################
# ITK
#####################
if ( NOT ANGIOTK_USE_SYSTEM_ITK )
  ExternalProject_Add(AngioTk_ExternalPackages_ITK    # Name for custom target
    PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/ITK"
    GIT_REPOSITORY "https://github.com/InsightSoftwareConsortium/ITK.git"
    GIT_TAG v4.11.1
    CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS} -DBUILD_EXAMPLES=OFF -DBUILD_TESTING=OFF -DITK_USE_REVIEW=ON -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/ITK
    #UPDATE_DISCONNECTED 1
    UPDATE_COMMAND ""
    )
  if (0)
    install(SCRIPT "${AngioTk_SOURCE_DIR}/CMake/AngioTkInstallExternalPackages.cmake")
  endif()
  set( ITK_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/ITK/lib/cmake/ITK-4.11)

  list(APPEND ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_ITK)
  list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_ITK)
endif()


#####################
# VMTK
#####################
if ( NOT ANGIOTK_USE_SYSTEM_VMTK )
  set( VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS})
  list( APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DVMTK_USE_SUPERBUILD=OFF)
  list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_ITK=ON -DITK_DIR=${ITK_DIR} )
  if ( FEELPP_HAS_PARAVIEW )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_VTK=ON -DParaView_DIR=${FEELPP_ParaView_DIR} -DVTK_DIR=${FEELPP_ParaView_DIR} )
  elseif( FEELPP_HAS_VTK )
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DUSE_SYSTEM_VTK=ON -DVTK_DIR=${FEELPP_VTK_DIR} )
  else()
  endif()
  IF (APPLE)
    list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_MACOSX_RPATH=1)
  endif()
  list(APPEND VMTK_SUPERBUILD_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/ExternalPackages/install/VMTK)

  ExternalProject_Add(AngioTk_ExternalPackages_VMTK    # Name for custom target
    PREFIX "${CMAKE_BINARY_DIR}/ExternalPackages/build/VMTK"
    GIT_REPOSITORY "https://github.com/vmtk/vmtk.git"
    GIT_TAG v1.3.2
    CMAKE_ARGS ${VMTK_SUPERBUILD_CMAKE_ARGS}
    #UPDATE_DISCONNECTED 1
    UPDATE_COMMAND ""
    DEPENDS ${ANGIOTK_VMTK_EXTERNALPROJECT_DEPENDS}
    )

  set(VMTK_DIR ${AngioTk_BINARY_DIR}/ExternalPackages/install/VMTK)
  list(APPEND ANGIOTK_EXTERNALPROJECT_DEPENDS AngioTk_ExternalPackages_VMTK)

endif()


#####################
# AngioTk
#####################


set( ANGIOTK_SUPERBUILD_CMAKE_ARGS
  -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER} -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}  -DBUILD_SHARED_LIBS=${BUILD_SHARED_LIBS}
  -DBUILD_EXAMPLES=${BUILD_EXAMPLES} -DBUILD_DOXYGEN=${BUILD_DOXYGEN}
  -DANGIOTK_USE_SUPERBUILD=OFF
  -DFEELPP_DIR=$ENV{FEELPP_DIR}
  -DANGIOTK_USE_SYSTEM_ITK=ON -DITK_DIR=${ITK_DIR}
  -DANGIOTK_USE_SYSTEM_VMTK=ON -DVMTK_DIR=${VMTK_DIR}
  -DCMAKE_INSTALL_PREFIX=${AngioTk_BINARY_DIR}/AngioTkSuperBuild/install
  )
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
