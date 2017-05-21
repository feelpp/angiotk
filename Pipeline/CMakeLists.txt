
set(ANGIOTK_EXECUTABLE_DIR ${AngioTk_BINARY_DIR}/bin)
configure_file(runAngioTKPipeline.py.in ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py DESTINATION bin)

configure_file(AngioTkInstallPipelineScript.cmake.in ${CMAKE_CURRENT_BINARY_DIR}/AngioTkInstallPipelineScript.cmake @ONLY)
INSTALL(SCRIPT  ${CMAKE_CURRENT_BINARY_DIR}/AngioTkInstallPipelineScript.cmake)

#configure_file(setupAngioTKEnvironment.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh DESTINATION bin)
