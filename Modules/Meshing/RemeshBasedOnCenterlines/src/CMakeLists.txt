
angiotk_add_component( meshing_centerlines centerlines_main.cpp )
target_link_libraries( meshing_centerlines ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})

angiotk_add_component( meshing_remeshstl remeshstl_main.cpp )
target_link_libraries( meshing_remeshstl ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})

angiotk_add_component( meshing_volume volume_main.cpp )
target_link_libraries( meshing_volume ${FEELPP_LIBRARY} ${FEELPP_APP_LINK_LIBRARIES} ${FEELPP_LIBRARIES})
