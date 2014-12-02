# this script finds the modules in Modules/ folder 
# and include the components for each enabled module

macro( AngioTkModuleEnablement _modulePath )
  SET( AngioTk_CURRENT_MODULE_PATH ${_modulePath})
  get_filename_component( moduleName ${_modulePath} NAME)
  SET( AngioTk_CURRENT_MODULE_NAME ${moduleName})
  SET( BUILD_MODULE_${moduleName} FALSE CACHE BOOL "Build module ${moduleName}")
 
  IF( ${BUILD_MODULE_${moduleName}} )
    angiotk_init_module_python(${moduleName}) 
    file( GLOB components ${_modulePath}/* ) 
    foreach( component ${components} )
      if( IS_DIRECTORY ${component} )
        ADD_SUBDIRECTORY( ${component} )  
      endif()
    endforeach()
  ENDIF()
endmacro() 


file( GLOB modules ${AngioTk_SOURCE_DIR}/Modules/* ) 
foreach( module ${modules} )
  if( IS_DIRECTORY ${module} )
    AngioTkModuleEnablement( ${module} )
  endif()
endforeach()
