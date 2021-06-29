########################################
# Set up how to compile the source files
########################################

# Add the source files
SET(IVEGen_src ${SRCIVEGen}/io.f90
               ${SRCIVEGen}/main.f90
               ${SRCIVEGen}/possu.f90
               ${SRCIVEGen}/rng.f90
               ${SRCIVEGen}/T_matrix.f90
               ${SRCIVEGen}/common.f90
               ${SRCIVEGen}/${RADIUSGEN}
               ${SRCIVEGen}/sfunctions.f90
               ${SRCIVEGen}/translations.f90
               ${SRCIVEGen}/geometry_x.f90
               ${SRCIVEGen}/octtree.f90
               ${SRCIVEGen}/solver.f90
               ${SRCIVEGen}/mie.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${IVEGenEXE} ${IVEGen_src})

#####################################################
# Add the needed libraries and special compiler flags
#####################################################

# Uncomment if you have parallization
include_directories(${HDF5_INCLUDE_DIR})
target_link_libraries(${IVEGenEXE} ${HDF5_Fortran_LIBRARIES})


# Uncomment if you have parallization
SET_TARGET_PROPERTIES(${IVEGenEXE} PROPERTIES
                      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(${IVEGenEXE} ${MPI_Fortran_LIBRARIES})

# Uncomment if you need to link to BLAS and LAPACK
tARGET_LINK_LIBRARIES(${IVEGenEXE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${CMAKE_THREAD_LIBS_INIT})
