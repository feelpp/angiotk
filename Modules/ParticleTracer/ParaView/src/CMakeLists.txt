
angiotk_add_component(angiotk_particletracer_pv particletracer_pv.cpp)
target_link_libraries(angiotk_particletracer_pv ${FEELPP_LIBRARIES})
#target_link_libraries(particletracer_pv ${PARAVIEW_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARY} ${Boost_LIBRARIES})
#install(TARGETS particletracer_pv RUNTIME DESTINATION Modules/ParticleTracer/bin)

add_library(angiotk_particletracer_H5ParticleReader SHARED H5ParticleReader.cpp)
target_link_libraries(angiotk_particletracer_H5ParticleReader ${FEELPP_LIBRARIES})
#target_link_libraries(H5ParticleReader ${PARAVIEW_LIBRARIES} ${VTK_LIBRARIES} ${HDF5_LIBRARY} ${Boost_LIBRARIES})
install(TARGETS angiotk_particletracer_H5ParticleReader LIBRARY DESTINATION lib COMPONENT Libs)
