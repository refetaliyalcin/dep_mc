#!bin/bash

cp makefiles/Makefile_rng src/dsfmt/Makefile
cp makefiles/Makefile_stmm src/stmm/Makefile
cp makefiles/Makefile_root Makefile

mkdir src/stmm -p
cp IVEGen/src/IVEGen/src/common.f90         src/stmm/
cp IVEGen/src/IVEGen/src/sfunctions.f90     src/stmm/
cp IVEGen/src/IVEGen/src/translations.f90   src/stmm/
cp IVEGen/src/IVEGen/src/mie.f90            src/stmm/
cp IVEGen/src/IVEGen/src/io.f90             src/stmm/
cp IVEGen/src/IVEGen/src/BHMIE.f90          src/stmm/
cp src/fullf/translations_extension.f90    src/stmm/
cp src/fullf/mie_extension.f90             src/stmm/
cp src/fullf/T_matrix.f90                  src/stmm/
cp src/error_handler.f90                    src/stmm
cp src/constants.f90                        src/stmm





##Preparations for R2T2

#Variable 1 = find me, Variable 2 = replace me, #Variable 3 = remove substring before the replace
replace_text () {
    var=$(grep -A1 "^$1$" src/stmm/Makefile | tail -n 1 | sed "s/$3*//g")
    #echo $var
    sed -i -e "s@$2@$var@g" Makefile
}



SRC_LOC="src"
STMM_FOLDER="stmm"
sed -i -e 's/#Makefile for RT-Engine/#Makefile for R2T2/g' Makefile
#Extract files from STMM folder
replace_text "#SCRIPT_FIND_ME_00520" "#REPLACE_ME_00520" "OBJECTS ="
replace_text "#I_REQUIRE_THESE_FLAGS_911" "#REPLACE_ME_911" "LIBS ="
sed -i 's|#REPLACE_DEFAULT_BUID|sphereR2T2|g' Makefile


