########################################
# Set up how to compile the source files
########################################

# Make sure the build type is uppercase
STRING(TOUPPER "${GEOMETRY}" BT)

IF(BT STREQUAL "SPHERE")
    MESSAGE(STATUS "SPHERE GEOMETRY" ...)
    SET (GEOMETRYFILE "geom_sphere.f90")
ELSEIF(BT STREQUAL "PLANE")
    MESSAGE(STATUS "PLANE GEOMETRY" ...)
    SET (GEOMETRYFILE "geom_plane.f90")
ELSE()
    MESSAGE(STATUS "SPHERE GEOMETRY" ...)
    SET (GEOMETRYFILE "geom_sphere.f90")
ENDIF(BT STREQUAL "SPHERE")




# Add the source files
SET(R2T2_src ${SRCR2T2}/absorption_module.f90
               ${SRCR2T2}/algorithm.f90
               ${SRCR2T2}/constants.f90
               ${SRCR2T2}/definitions.f90
               ${SRCR2T2}/error_handler.f90
               ${SRCR2T2}/${GEOMETRYFILE}
               ${SRCR2T2}/grid_holder.f90
               ${SRCR2T2}/hashtable.f90
               ${SRCR2T2}/inputparser.f90
               ${SRCR2T2}/input_reader.f90
               ${SRCR2T2}/main.f90
               ${SRCR2T2}/math_routines.f90
               ${SRCR2T2}/misc.f90
               ${SRCR2T2}/mpi_tools.f90
               ${SRCR2T2}/output_writer.f90
               ${SRCR2T2}/preparedata.f90
               ${SRCR2T2}/scatterer_r2t2.f90
               ${SRCR2T2}/sfmtrng.f90
)

# Define the executable in terms of the source files
ADD_EXECUTABLE(${R2T2EXE} ${R2T2_src})



TARGET_LINK_LIBRARIES(${R2T2EXE} ${DSFMTLIB})
TARGET_LINK_LIBRARIES(${R2T2EXE} ${FULLFLIB})


#####################################################
# Add the needed libraries and special compiler flags
#####################################################


# Uncomment if you have parallization
SET_TARGET_PROPERTIES(${R2T2EXE} PROPERTIES
                      COMPILE_FLAGS "${MPI_Fortran_COMPILE_FLAGS}"
                      LINK_FLAGS "${MPI_Fortran_LINK_FLAGS}")
INCLUDE_DIRECTORIES(${MPI_Fortran_INCLUDE_PATH})
TARGET_LINK_LIBRARIES(${R2T2EXE} ${MPI_Fortran_LIBRARIES})

# Uncomment if you need to link to BLAS and LAPACK
tARGET_LINK_LIBRARIES(${R2T2EXE} ${BLAS_LIBRARIES}
                                ${LAPACK_LIBRARIES}
                                ${HDF5_Fortran_LIBRARIES})

