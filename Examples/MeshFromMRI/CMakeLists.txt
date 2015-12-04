
# Get all the examples subdirectories
# and add them to the cmake process
file(GLOB files "*")
foreach(fl ${files})
    if(IS_DIRECTORY ${fl} AND ${fl} MATCHES "Base2_FOA42")
        add_subdirectory(${fl})
    endif()
endforeach()
