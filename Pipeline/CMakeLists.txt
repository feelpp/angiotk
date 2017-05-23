
set(ANGIOTK_EXECUTABLE_DIR ${AngioTk_BINARY_DIR}/bin)
configure_file(angiotk_pipeline_run.py.in ${CMAKE_CURRENT_BINARY_DIR}/angiotk_pipeline_run.py)

install(PROGRAMS  ${CMAKE_CURRENT_BINARY_DIR}/angiotk_pipeline_run.py
  DESTINATION bin )  
