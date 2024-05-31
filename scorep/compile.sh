#!/bin/bash

if [ $# -lt 1 ]
then
    printf "\n\t Usage: ./compile.sh <compiler (e.g. gcc, nvcc, ecc.)> \n\n"
    exit 1
else
    CC=$1
    if ! command -v ${CC} &>/dev/null
    then
	printf "\n\t ${CC} not found... aborting... \n"
	exit 2
    fi

    if ! command -v scorep &>/dev/null
    then
	printf "\n\t scorep not found... aborting... \n"
	exit 3
    fi
fi

printf "\n\n\t Using the ${CC} compiler \n\n"

cd ${WORKDIR}/src

for EXE in ${EXEC[@]}
do
    DEBUG_FLAG=$(echo ${EXE} | grep "Debug")
    OMP_FLAG=$(echo ${EXE} | grep "OMP")

    if [ ${#DEBUG_FLAG} -eq 0 ]
    then
	DEBUG_SWITCH=NO
    else
	DEBUG_SWITCH=YES
    fi

    if [ ${#OMP_FLAG} -eq 0 ]
    then
	OMP_SWITCH=NO
    else
	OMP_SWITCH=YES
    fi

    # printf "\n\t Compiling ${EXE}... \n"
    # COMPILE=compile_${EXE}.txt
    # make clean && make EXEC=${EXE} COMPILER=${CC} DEBUG=${DEBUG} OMP=${OMP} SCOREP=NO |& tee ${COMPILE}
    # file ${EXE} |& tee -a ${COMPILE}
    # ldd ${EXE} |& tee -a ${COMPILE}
    # mv ${EXE} ${COMPILE} ../example

    # Scorep compilation
    EXE_SCOREP=${EXE}Scorep
    printf "\n\t Compiling ${EXE}... \n"
    COMPILE=compile_${EXE}.txt
    make clean && make EXEC=${EXE_SCOREP} COMPILER=${CC} DEBUG=${DEBUG_SWITCH} OMP=${OMP_SWITCH} SYSTYPE=${SYSTEM} SCOREP=YES |& tee ${COMPILE}
    file ${EXE_SCOREP} |& tee -a ${COMPILE}
    ldd ${EXE_SCOREP} |& tee -a ${COMPILE}
    mv ${EXE_SCOREP} ${COMPILE} ${OUT_DIR}

done

make clean
cd ${WORKDIR}

printf "\n\t Compilation done \n"
