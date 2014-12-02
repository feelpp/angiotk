macro( angiotk_init_module_python _moduleName )
  FILE( WRITE ${AngioTk_BINARY_DIR}/${_moduleName}.py
        "import os\n\n"
      )
endmacro()

macro( angiotk_add_component _target _srcs )
  message( STATUS ${AngioTk_CURRENT_MODULE_NAME}.${_target})
  add_executable( ${_target} ${_srcs} )
  FILE( APPEND ${AngioTk_BINARY_DIR}/${AngioTk_CURRENT_MODULE_NAME}.py
        "def ${_target}(input, output):\n"
        "  os.system(\"${AngioTk_BINARY_DIR}/bin/${_target} \"+input+\" \"+output)\n\n"
      )
endmacro()

  
