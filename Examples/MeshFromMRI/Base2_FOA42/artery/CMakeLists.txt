
# Compute the relative path from top source dir
#message(${CMAKE_SOURCE_DIR})
#message(${CMAKE_CURRENT_SOURCE_DIR})
string(REPLACE "${CMAKE_SOURCE_DIR}/" "" RELATIVE_PATH_FROM_TOP_SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR})
#message(${RELATIVE_PATH_FROM_TOP_SOURCE_DIR})

# Configure and install the script for executing the pipeline
configure_file(runpipeline.sh.in ${CMAKE_CURRENT_BINARY_DIR}/runpipeline.sh @ONLY)
install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runpipeline.sh
    DESTINATION share/AngioTK/${RELATIVE_PATH_FROM_TOP_SOURCE_DIR})

# Install the cfg files
file(GLOB cfgfiles "*.cfg")
install(FILES ${cfgfiles}
    DESTINATION share/AngioTK/${RELATIVE_PATH_FROM_TOP_SOURCE_DIR})

# Copy the data directory
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
    DESTINATION share/AngioTK/${RELATIVE_PATH_FROM_TOP_SOURCE_DIR})

