# Get AngioTK "version number": the last commit ID (abbreviated)
# Store it as CMAKE_ANGIOTK_VERSION
execute_process(
    COMMAND git rev-parse --verify -q --short HEAD
    WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
    OUTPUT_VARIABLE CMAKE_ANGIOTK_VERSION
    ERROR_VARIABLE GIT_COMMIT_error
    RESULT_VARIABLE GIT_COMMIT_result)
  if ( NOT ${GIT_COMMIT_result} EQUAL 0 )
    MESSAGE(SEND_ERROR "Command \"git rev-parse --verify -q --short git-svn\" failed with output:\n${GIT_COMMIT_error}")
  else()
    STRING(STRIP "${CMAKE_ANGIOTK_VERSION}" CMAKE_ANGIOTK_VERSION )
  endif()

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

configure_file(screenshot.py.in ${CMAKE_CURRENT_BINARY_DIR}/screenshot.py)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/screenshot.py DESTINATION bin)
