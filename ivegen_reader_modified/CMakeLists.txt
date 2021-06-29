# CMake project file for Incoherent Volume Element Generator (STMM)
#
# Original file by SethMMorton
# https://github.com/SethMMorton/cmake_fortran_template
##################################################
# Define the project and the depencies that it has
##################################################

CMAKE_MINIMUM_REQUIRED(VERSION 3.0.2)
PROJECT(IVEGen Fortran)

# Set the IVEGen version
SET(VERSION 1.0)

# Add our local modules to the module path
SET(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/Modules/")

# If argument USE_POWERLAW_DISTRIBUTION is given, use a
# power law to generate radiuses
SET(USE_DISTRIBUTION "normal" CACHE STRING "Some user-specified option")
IF(USE_DISTRIBUTION STREQUAL "powerlaw")
    MESSAGE("The radiuses of the spheres will obey a power law distribution.")
    SET(RADIUSGEN radiusgenpower.f90)
ELSEIF(USE_DISTRIBUTION STREQUAL "normal")
    MESSAGE("The radiuses of the spheres will obey a normal distribution.")
    SET(RADIUSGEN radiusgennormal.f90)
ELSE()
    MESSAGE(FATAL_ERROR "Invalid distribution")
ENDIF()



# Uncomment if it is required that Fortran 90 is supported
IF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)
    MESSAGE(FATAL_ERROR "Fortran compiler does not support F90")
ENDIF(NOT CMAKE_Fortran_COMPILER_SUPPORTS_F90)

# Set some options the user may choose
# Uncomment the below if you want the user to choose a parallelization library
IF (NOT MPI_Fortran_FOUND)
    FIND_PACKAGE (MPI REQUIRED)
    SET(USE_MPI)
ENDIF (NOT MPI_Fortran_FOUND)


# Setup the LAPACK libraries.  This also takes care of peculiarities, such as
# the fact the searching for MKL requires a C compiler, and that the results
# are not stored in the cache. 
ENABLE_LANGUAGE(C) # Some libraries need a C compiler to find

FIND_PACKAGE(LAPACK REQUIRED)
# Remember that LAPACK (and BLAS) was found.  For some reason the
# FindLAPACK routine doesn't place these into the CACHE.
SET(BLAS_FOUND TRUE CACHE INTERNAL "BLAS was found" FORCE)
SET(LAPACK_FOUND TRUE CACHE INTERNAL "LAPACK was found" FORCE)
SET(BLAS_LIBRARIES ${BLAS_LIBRARIES} CACHE INTERNAL "BLAS LIBS" FORCE)
SET(LAPACK_LIBRARIES ${LAPACK_LIBRARIES} CACHE INTERNAL "LAPACK LIBS" FORCE)
SET(BLA_STATIC TRUE)
#Require HDF5
FIND_PACKAGE(HDF5 COMPONENTS Fortran REQUIRED)


# This INCLUDE statement executes code that sets the compile flags for DEBUG,
# RELEASE, and TESTING.  You should  review this file and make sure the flags 
# are to your liking.
INCLUDE(${CMAKE_MODULE_PATH}/SetFortranFlags.cmake) 


# There is an error in CMAKE with this flag for pgf90.  Unset it
GET_FILENAME_COMPONENT(FCNAME ${CMAKE_Fortran_COMPILER} NAME)
IF(FCNAME STREQUAL "pgf90")
    UNSET(CMAKE_SHARED_LIBRARY_LINK_Fortran_FLAGS)
ENDIF(FCNAME STREQUAL "pgf90")





############################################################
# Define the actual files and folders that make up the build
############################################################

# Define the executable name
SET(IVEGenEXE IVEGen)

# Define some directories
SET(SRC ${CMAKE_SOURCE_DIR}/src)
SET(LIB ${CMAKE_SOURCE_DIR}/lib)
SET(BIN ${CMAKE_SOURCE_DIR}/bin)
SET(SRCIVEGen ${SRC})

# Have the .mod files placed in the lib folder
SET(CMAKE_Fortran_MODULE_DIRECTORY ${LIB})

# The source for the IVEGen binary and have it placed in the bin folder
ADD_SUBDIRECTORY(${SRCIVEGen} ${BIN})

# Add a distclean target to the Makefile
ADD_CUSTOM_TARGET(distclean 
    COMMAND ${CMAKE_COMMAND} -P ${CMAKE_SOURCE_DIR}/distclean.cmake
)


