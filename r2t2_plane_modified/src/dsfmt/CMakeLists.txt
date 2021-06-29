ENABLE_LANGUAGE(C) 
SET(DSFMT_src ${SRCDSFMT}/dSFMT_interface.F90
              ${SRCDSFMT}/dsfmt/dSFMT.c
              ${SRCDSFMT}/dsfmt/dSFMT-jump.c
              ${SRCDSFMT}/dsfmt/dSFMT_utils.c
)

##
set(CMAKE_Fortran_FLAGS_RELEASE "")
set(CMAKE_C_FLAGS_RELEASE "")

#Set DSFMT_MEXP
set(CMAKE_Fortran_FLAGS "-DDSFMT_MEXP=1279")
set(CMAKE_C_FLAGS  "-DDSFMT_MEXP=1279")

##Flags for fortran compiler (GCC)
if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -cpp ")
endif()

if("${CMAKE_Fortran_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_Fortran_FLAGS  "${CMAKE_Fortran_FLAGS} -cpp ")
endif()


##Flags for C compiler (GCC)
if("${CMAKE_C_COMPILER_ID}" STREQUAL "GNU")
    #add_compile_options("-finline-functions" "-fomit-frame-pointer" "-DNDEBUG" "-fno-strict-aliasing" "--param max-inline-insns-single=1800" "-Wmissing-prototypes" "-Wall"  "-std=c99" "--param inline-unit-growth=500" "--param large-function-growth=900" "-msse2 -DHAVE_SSE2")
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}  -O3 -finline-functions -fomit-frame-pointer -DNDEBUG -fno-strict-aliasing --param max-inline-insns-single=1800 -Wall -std=c99 --param inline-unit-growth=500 --param large-function-growth=900 -msse2 -DHAVE_SSE2 ")
endif()

##Flags for C compiler (GCC)
if("${CMAKE_C_COMPILER_ID}" STREQUAL "Intel")
    set(CMAKE_C_FLAGS  "${CMAKE_C_FLAGS}  -O3 -finline-functions -fomit-frame-pointer -DNDEBUG  -fno-strict-aliasing -std=c99")
endif()


ADD_LIBRARY(${DSFMTLIB} STATIC ${DSFMT_src})



