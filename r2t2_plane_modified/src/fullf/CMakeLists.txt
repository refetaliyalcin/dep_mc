########################################
# Set up how to compile the source files
########################################





include(ExternalProject)
ExternalProject_Add(IVEGen
    GIT_REPOSITORY git@bitbucket.org:planetarysystemresearch/ivegen_pub.git
    PREFIX ${CMAKE_SOURCE_DIR}/IVEGen/
    DOWNLOAD_DIR ${CMAKE_SOURCE_DIR}/IVEGen/
    CONFIGURE_COMMAND ""
    BUILD_COMMAND     ""
    INSTALL_COMMAND   ""
    TEST_COMMAND      ""
)

SET(IVEGENDIR ${CMAKE_SOURCE_DIR}/IVEGen/src/IVEGen)


SET(COMMON ${IVEGENDIR}/src/common.f90)

SET(FULLF_src   ${COMMON}
                ${CMAKE_SOURCE_DIR}/IVEGen/src/IVEGen/src/sfunctions.f90
                ${CMAKE_SOURCE_DIR}/IVEGen/src/IVEGen/src/translations.f90
                ${CMAKE_SOURCE_DIR}/IVEGen/src/IVEGen/src/mie.f90
                ${CMAKE_SOURCE_DIR}/IVEGen/src/IVEGen/src/io.f90
                ${SRCFULLF}/translations_extension.f90
                ${SRCFULLF}/mie_extension.f90
                ${SRCFULLF}/T_matrix.f90
)


if (NOT EXISTS "${COMMON}")
    execute_process(COMMAND "git" "clone" "git@bitbucket.org:planetarysystemresearch/ivegen_pub.git" "${IVEGENDIR}")
endif()

ADD_LIBRARY(${FULLFLIB} STATIC ${FULLF_src})

include_directories(${HDF5_INCLUDE_DIR})
target_link_libraries(${FULLFLIB} ${HDF5_Fortran_LIBRARIES})

add_dependencies(${FULLFLIB} IVEGen)
