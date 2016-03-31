
add_executable(test_reader test_reader.cpp)
target_link_libraries(test_reader ${PARAVIEW_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARY} ${Boost_LIBRARIES} H5ParticleReader)
install(TARGETS test_reader RUNTIME DESTINATION Modules/ParticleTracer/bin)
