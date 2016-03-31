
add_executable(particletracer_pv particletracer_pv.cpp)
target_link_libraries(particletracer_pv ${PARAVIEW_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARY} ${Boost_LIBRARIES})
install(TARGETS particletracer_pv RUNTIME DESTINATION Modules/ParticleTracer/bin)

add_library(H5ParticleReader H5ParticleReader.cpp)
target_link_libraries(H5ParticleReader ${PARAVIEW_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARY} ${Boost_LIBRARIES})
install(TARGETS H5ParticleReader LIBRARY DESTINATION Modules/ParticleTracer/lib)
