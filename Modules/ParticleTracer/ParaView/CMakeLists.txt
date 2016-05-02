# Module ParticleTracer
# Label volume to Mesh for CFD

project(PVParticleTracer)

message(STATUS "Configuring ParticleTracer (ParaView) Module")

# reset all include directories to avoid errors with different VTK versions
# We want to use the version provided by ParaView
set_property(DIRECTORY . PROPERTY INCLUDE_DIRECTORIES "")

# Don't use the basic FindHDF5.cmake package
# It often ends up mixing custom installations and system installations
#FIND_PACKAGE(HDF5 REQUIRED)

find_path(HDF5_INCLUDE_DIR hdf5.h
    PATHS $ENV{HDF5_DIR}/include
    NO_SYSTEM_PATH
    )

find_library(HDF5_LIBRARY hdf5
    PATHS $ENV{HDF5_DIR}/lib
    NO_SYSTEM_PATH
    )

if(HDF5_INCLUDE_DIR AND HDF5_LIBRARY)
    include_directories(${HDF5_INCLUDE_DIR})
else()
    message(FATAL_ERROR "Couldn't find HDF5 libraries. Please set HDF5_DIR environment variables.")
endif()

find_package(Boost REQUIRED COMPONENTS program_options)

if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIR} ${BOOST_INCLUDE_PATH})
endif()

FIND_PACKAGE(ParaView REQUIRED
    COMPONENTS vtkParallelMPI vtkPVCatalyst vtkPVPythonCatalyst
    PATHS $ENV{PARAVIEW_DIR} ${MACHINE_PARAVIEW_DIR})

INCLUDE(${PARAVIEW_USE_FILE})

include_directories(include)

add_subdirectory(src)
add_subdirectory(test)

