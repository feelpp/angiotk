
set(ANGIOTK_EXECUTABLE_DIR ${AngioTk_BINARY_DIR}/bin)
configure_file(angiotk_pipeline_run.py.in ${CMAKE_CURRENT_BINARY_DIR}/angiotk_pipeline_run.py)

configure_file(AngioTkInstallPipelineScript.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/AngioTkInstallPipelineScript.cmake @ONLY)
INSTALL(SCRIPT  ${CMAKE_CURRENT_BINARY_DIR}/AngioTkInstallPipelineScript.cmake)

#configure_file(setupAngioTKEnvironment.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh DESTINATION bin)
