macro( angiotk_add_component _target _srcs )
  message( STATUS ${AngioTk_CURRENT_MODULE_NAME}.${_target})
  add_executable( ${_target} ${_srcs} )
endmacro()

  
