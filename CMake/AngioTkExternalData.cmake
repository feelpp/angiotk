get_filename_component(_AngioTkExternalData_DIR "${CMAKE_CURRENT_LIST_FILE}" PATH)
include(${_AngioTkExternalData_DIR}/ExternalData.cmake)

set(ExternalData_BINARY_ROOT ${CMAKE_BINARY_DIR}/ExternalData)

set(ExternalData_URL_TEMPLATES "" CACHE STRING
  "Additional URL templates for the ExternalData CMake script to look for testing data. E.g.
file:///var/bigharddrive/%(algo)/%(hash)")
mark_as_advanced(ExternalData_URL_TEMPLATES)
list(APPEND ExternalData_URL_TEMPLATES
  # Data published by MIDAS
  "http://vivabrain.u-strasbg.fr/midas/api/rest?method=midas.bitstream.download&checksum=%(hash)&algorithm=%(algo)"
  )

# Tell ExternalData commands to transform raw files to content links.
# TODO: Condition this feature on presence of our pre-commit hook.
set(ExternalData_LINK_CONTENT MD5)


#-----------------------------------------------------------------------------
# AngioTk wrapper for add_test that automatically sets the test's LABELS property
# to the value of its containing module.
#
function(angiotk_add_test)
  # Add tests with data in the ITKData group.
  ExternalData_add_test(angioTkData ${ARGN})
endfunction()

# Macro to create a test driver
macro(CreateTestDriver NAME LIBS KitTests)
  set( ADDITIONAL_SRC ${ARGN} )
  set(CMAKE_TESTDRIVER_BEFORE_TESTMAIN "#include \"itkTestDriverBeforeTest.inc\"")
  set(CMAKE_TESTDRIVER_AFTER_TESTMAIN "#include \"itkTestDriverAfterTest.inc\"")
  create_test_sourcelist(Tests ${NAME}TestDriver.cxx
    ${KitTests}
	EXTRA_INCLUDE itkTestDriverIncludeRequiredIOFactories.h
    FUNCTION  ProcessArgumentsAndRegisterRequiredFactories
    )
  add_executable(${NAME}TestDriver ${NAME}TestDriver.cxx ${Tests} ${ADDITIONAL_SRC})
  target_link_libraries(${NAME}TestDriver ${LIBS})
endmacro()
