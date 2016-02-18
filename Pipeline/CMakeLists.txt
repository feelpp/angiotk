
configure_file(runAngioTKPipeline.py.in ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py DESTINATION bin)
