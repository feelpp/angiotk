message(STATUS "[angiotk] adding missing Gmsh headers...")

# we include this directory : add some missing headers from Gmsh
find_program(WGET NAMES wget )
if ( WGET )
  add_custom_target(angiotk_gmsh_get_files ALL)
  set(gmsh_gitlab_url "http://gitlab.onelab.info/gmsh/gmsh/raw/gmsh_${GMSH_MAJOR_VERSION}_${GMSH_MINOR_VERSION}_${GMSH_PATCH_VERSION}")
  set( GMSH_VERSION_FULL ${GMSH_MAJOR_VERSION}.${GMSH_MINOR_VERSION}.${GMSH_PATCH_VERSION})
  if ( GMSH_VERSION_FULL VERSION_LESS 2.16 )
    add_custom_command(
      TARGET angiotk_gmsh_get_files
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/Field.h -O Field.h 
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Numeric/BasisFactory.h -O BasisFactory.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Common/OpenFile.h -O OpenFile.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/BackgroundMesh.h -O BackgroundMesh.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Geo/Curvature.h -O Curvature.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/meshPartitionObjects.h -O meshPartitionObjects.h
      )
    # install(FILES
    #   ${CMAKE_CURRENT_BINARY_DIR}/BasisFactory.h ${CMAKE_CURRENT_BINARY_DIR}/Field.h ${CMAKE_CURRENT_BINARY_DIR}/OpenFile.h
    #   ${CMAKE_CURRENT_BINARY_DIR}/BackgroundMesh.h ${CMAKE_CURRENT_BINARY_DIR}/Curvature.h ${CMAKE_CURRENT_BINARY_DIR}/meshPartitionObjects.h
    #   DESTINATION include/feelpp/gmsh COMPONENT Devel )
  else()
    add_custom_command(
      TARGET angiotk_gmsh_get_files
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/Field.h -O Field.h 
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Numeric/BasisFactory.h -O BasisFactory.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Common/OpenFile.h -O OpenFile.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/BackgroundMesh.h -O BackgroundMesh.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Geo/Curvature.h -O Curvature.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/meshPartitionObjects.h -O meshPartitionObjects.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/BackgroundMeshTools.h -O BackgroundMeshTools.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/qualityMeasuresJacobian.h -O qualityMeasuresJacobian.h
      COMMAND ${WGET} -N -q --timeout=5 ${gmsh_gitlab_url}/Mesh/meshGRegion.h -O meshGRegion.h
      )
    # install(FILES
    #   ${CMAKE_CURRENT_BINARY_DIR}/BasisFactory.h ${CMAKE_CURRENT_BINARY_DIR}/Field.h ${CMAKE_CURRENT_BINARY_DIR}/OpenFile.h
    #   ${CMAKE_CURRENT_BINARY_DIR}/BackgroundMesh.h ${CMAKE_CURRENT_BINARY_DIR}/Curvature.h ${CMAKE_CURRENT_BINARY_DIR}/meshPartitionObjects.h
    #   ${CMAKE_CURRENT_BINARY_DIR}/BackgroundMeshTools.h ${CMAKE_CURRENT_BINARY_DIR}/qualityMeasuresJacobian.h ${CMAKE_CURRENT_BINARY_DIR}/meshGRegion.h
    #   DESTINATION include/feelpp/gmsh COMPONENT Devel )
  endif()
endif()

