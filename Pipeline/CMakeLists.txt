
configure_file(runAngioTKPipeline.py.in ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runAngioTKPipeline.py DESTINATION bin)

configure_file(setupAngioTKEnvironment.sh.in ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/setupAngioTKEnvironment.sh DESTINATION bin)

configure_file(master.py.in ${CMAKE_CURRENT_BINARY_DIR}/master.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/master.py DESTINATION bin)

configure_file(testAngioTKPipeline.py.in ${CMAKE_CURRENT_BINARY_DIR}/testAngioTKPipeline.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/testAngioTKPipeline.py DESTINATION bin)

configure_file(createAsciidocDB.py.in ${CMAKE_CURRENT_BINARY_DIR}/createAsciidocDB.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/createAsciidocDB.py DESTINATION bin)

configure_file(createAsciidocSummary.py.in ${CMAKE_CURRENT_BINARY_DIR}/createAsciidocSummary.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/createAsciidocSummary.py DESTINATION bin)

configure_file(screenShot.py.in ${CMAKE_CURRENT_BINARY_DIR}/screenShot.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/screenShot.py DESTINATION bin)

configure_file(screenshot.py.in ${CMAKE_CURRENT_BINARY_DIR}/screenshot.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/screenshot.py DESTINATION bin)
