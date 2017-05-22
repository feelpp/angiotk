

set(ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR ${PostProcessing_BINARY_DIR}/test)
set(ExternalData_SOURCE_ROOT ${PostProcessing_SOURCE_DIR})
set(ExternalData_BINARY_ROOT ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR})
set(EXTERNALDATA_ALL_OUTPUT_FILES)
ExternalData_expand_arguments( angioTkData EXTERNALDATA_OUTPUT_FILES DATA{${PostProcessing_SOURCE_DIR}/simpleTbifurcationFSI.msh} )
list(APPEND EXTERNALDATA_ALL_OUTPUT_FILES ${EXTERNALDATA_OUTPUT_FILES})
unset(ExternalData_BINARY_ROOT)
unset(ExternalData_SOURCE_ROOT)

# install data downloaded (special treatment (data are not generated at this time and file are symbolic link)
install(CODE 
  "foreach( datafile ${EXTERNALDATA_ALL_OUTPUT_FILES})
     get_filename_component( _realfile \"\${datafile}\" REALPATH)
     get_filename_component( _name \"\${datafile}\" NAME)
     file(INSTALL \${_realfile} DESTINATION \${CMAKE_INSTALL_PREFIX}/share/AngioTK/Modules/Meshing/PostProcessing/Tests RENAME \${_name} )
   endforeach()"
)

add_dependencies(angiotk_meshing_extractsubmesh angioTkData)
add_dependencies(angiotk_meshing_meshpartitioner angioTkData)
add_dependencies(angiotk_meshing_postprocessing angioTkData)


set(ANGIOTK_MODULE_MESHING_TEST_CFG_FILES extractsubmesh.cfg meshpartitioner.cfg postprocessing.cfg)
set(ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL)
foreach(cfgfile ${ANGIOTK_MODULE_MESHING_TEST_CFG_FILES})
  configure_file( ${cfgfile} ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/${cfgfile} COPYONLY)
  list(APPEND ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/${cfgfile})
endforeach()
install(FILES ${ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL} DESTINATION share/AngioTK/Modules/Meshing/PostProcessing/Tests)



add_test(NAME test_angiotk_meshing_extractsubmesh COMMAND angiotk_meshing_extractsubmesh
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/extractsubmesh.cfg )

add_test(NAME test_angiotk_meshing_meshpartitioner COMMAND angiotk_meshing_meshpartitioner
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/meshpartitioner.cfg )

add_test(NAME test_angiotk_meshing_postprocessing COMMAND angiotk_meshing_pipeline_postprocessing 
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/postprocessing.cfg )
