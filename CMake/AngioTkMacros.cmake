# Given the absolute path to a module, this macro set the 
# AngioTk_CURRENT_MODULE_PATH/NAME variables and cache the 
# BUILD_MODULE_xxxx user option
macro( angiotk_module_enablement _modulePath )
  set( AngioTk_CURRENT_MODULE_PATH ${_modulePath})
  get_filename_component( moduleName ${_modulePath} NAME)
  set( AngioTk_CURRENT_MODULE_NAME ${moduleName})
  set( BUILD_MODULE_${moduleName} FALSE CACHE BOOL "Build module ${moduleName}")
endmacro()

# Initialisation of a python file containing definitions
# of components for a given module _moduleName
macro( angiotk_init_module_python _moduleName )
  file( WRITE ${AngioTk_EXECUTABLE_PATH}/${_moduleName}.py
        "# AngioTk components for ${_moduleName} module\n"
        "# This file has been generated automatically by CMake.\n"
        "# Any change will be overwritten on next Configure.\n"
        "import os\n\n"
      )
endmacro()

# Create a new AngioTk component:
#   - Add an executable _target from sources _srcs
#   - Define a function in the current module python file
# Note1: This should be called from inside a module (ie 
# ${AngioTk_CURRENT_MODULE_NAME} should be defined)
# Note2: the executable should take at least 2 args: 
#   - an input filename
#   - an output filename
# a third argument could optionnally be used to give additional
# component-specific arguments 
macro( angiotk_add_component _target _srcs )
  message( STATUS "Creating python file for ${AngioTk_CURRENT_MODULE_NAME}.${_target}")
  add_executable( ${_target} ${_srcs} )
  file( APPEND ${AngioTk_EXECUTABLE_PATH}/${AngioTk_CURRENT_MODULE_NAME}.py
        "def ${_target}(inputFile, outputFile, componentArgs=\"\"):\n"
        "  os.system(\"${AngioTk_EXECUTABLE_PATH}/${_target} \"
                     +inputFile+\" \"
                     +outputFile+\" \"
                     +componentArgs
                    )\n\n"
      )
endmacro()

  
