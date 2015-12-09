
# Configure and install the script for executing the pipeline
#configure_file(runpipeline.sh.in ${CMAKE_CURRENT_BINARY_DIR}/runpipeline.sh @ONLY)
#install(FILES ${CMAKE_CURRENT_BINARY_DIR}/runpipeline.sh
    #DESTINATION share/AngioTK/Examples/MeshFromMRI/Base2_GEO26/veinous)

# Install the cfg files
file(GLOB cfgfiles "*.cfg")
install(FILES ${cfgfiles}
    DESTINATION share/AngioTK/Examples/MeshFromMRI/Base2_GEO26/veinous)

# Copy the data directory
install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/data
    DESTINATION share/AngioTK/Examples/MeshFromMRI/Base2_GEO26/veinous)
