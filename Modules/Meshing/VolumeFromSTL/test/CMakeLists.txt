# Module Meshing
# Label volume to Mesh for CFD

set(ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR ${VolumeFromSTL_BINARY_DIR}/test)
set(ExternalData_SOURCE_ROOT ${VolumeFromSTL_SOURCE_DIR})
set(ExternalData_BINARY_ROOT ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR})
set(EXTERNALDATA_ALL_OUTPUT_FILES)
ExternalData_expand_arguments( angioTkData EXTERNALDATA_OUTPUT_FILES DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation.stl} )
list(APPEND EXTERNALDATA_ALL_OUTPUT_FILES ${EXTERNALDATA_OUTPUT_FILES})
ExternalData_expand_arguments( angioTkData EXTERNALDATA_OUTPUT_FILES DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation_remeshVMTKarea0.005.stl} )
list(APPEND EXTERNALDATA_ALL_OUTPUT_FILES ${EXTERNALDATA_OUTPUT_FILES})
ExternalData_expand_arguments( angioTkData EXTERNALDATA_OUTPUT_FILES DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation_remeshVMTKarea0.005_centerlines.vtk} )
list(APPEND EXTERNALDATA_ALL_OUTPUT_FILES ${EXTERNALDATA_OUTPUT_FILES})
ExternalData_expand_arguments( angioTkData EXTERNALDATA_OUTPUT_FILES DATA{${VolumeFromSTL_SOURCE_DIR}/simpleTbifurcation.inletoutlet.desc} )
list(APPEND EXTERNALDATA_ALL_OUTPUT_FILES ${EXTERNALDATA_OUTPUT_FILES})
unset(ExternalData_BINARY_ROOT)
unset(ExternalData_SOURCE_ROOT)

# install data downloaded (special treatment (data are not generated at this time and file are symbolic link)
install(CODE 
  "foreach( datafile ${EXTERNALDATA_ALL_OUTPUT_FILES})
     get_filename_component( _realfile \"\${datafile}\" REALPATH)
     get_filename_component( _name \"\${datafile}\" NAME)
     file(INSTALL \${_realfile} DESTINATION \${CMAKE_INSTALL_PREFIX}/share/AngioTK/Modules/Meshing/VolumeFromSTL/Tests RENAME \${_name} )
   endforeach()"
)

add_dependencies( angiotk_meshing_pipeline_volumefromstl angioTkData)

set(ANGIOTK_MODULE_MESHING_TEST_CFG_FILES centerlines.cfg remeshstlvmtk.cfg remeshstlgmsh.cfg volumefromstlandcenterlines.cfg volumefromstl.cfg)
set(ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL)
foreach(cfgfile ${ANGIOTK_MODULE_MESHING_TEST_CFG_FILES})
  configure_file( ${cfgfile} ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/${cfgfile} COPYONLY)
  list(APPEND ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL ${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/${cfgfile})
endforeach()
install(FILES ${ANGIOTK_MODULE_MESHING_TEST_CFG_FILES_TO_INSTALL} DESTINATION share/AngioTK/Modules/Meshing/VolumeFromSTL/Tests)


add_test(NAME test_angiotk_meshing_centerlines COMMAND angiotk_meshing_centerlines
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/centerlines.cfg )

add_test(NAME test_angiotk_meshing_remeshstlvmtk COMMAND angiotk_meshing_remeshstl
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/remeshstlvmtk.cfg )

add_test(NAME test_angiotk_meshing_remeshstlgmsh COMMAND angiotk_meshing_remeshstl
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/remeshstlgmsh.cfg )

add_test(NAME test_angiotk_meshing_volumefromstlandcenterlines COMMAND angiotk_meshing_volumefromstlandcenterlines
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/volumefromstlandcenterlines.cfg )

add_test(NAME test_angiotk_meshing_volumefromstl COMMAND angiotk_meshing_pipeline_volumefromstl
  --config-file=${ANGIOTK_MODULE_MESHING_TEST_BINARY_DIR}/volumefromstl.cfg )
