#!/bin/bash

if [ $# -lt 2 ]
then
    printf "\n\t Usage: ./compile.sh <C compiler> <CPP compiler> \n\n"
    exit 1
else
    CC=$1
    if ! command -v ${CC} &>/dev/null
    then
	printf "\n\t ${CC} not found... aborting... \n"
	exit 2
    fi

    if [ ${ENERGY} -eq 1 ]
    then
	CPP=$2
	if ! command -v ${CPP} &>/dev/null
	then
	    printf "\n\t ${CPP} not found... aborting... \n"
	    exit 3
	fi
    fi

    if [ ${PROFILING} -eq 1 ]
    then
	if ! command -v scorep &>/dev/null
	then
	    printf "\n\t scorep not found... aborting... \n"
	    exit 4
	fi
    fi
fi

printf "\n\n\t Using the ${CC} compiler \n\n"

cd ${WORKDIR}/src

for EXE in ${EXEC[@]}
do
    DEBUG_FLAG=$(echo ${EXE}  | grep "Debug")
    OMP_FLAG=$(echo ${EXE}    | grep "OMP")
    GPU_FLAG=$(echo ${EXE}    | grep "GPU")

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

    if [ ${#GPU_FLAG} -eq 0 ]
    then
	GPU_SWITCH=NO
    else
	GPU_SWITCH=YES
    fi
    
    if [ ${PROFILING} -eq 1 ] || [ ${TRACING} -eq 1 ]
    then
	# Scorep compilation
	EXE_=${EXE}_Scorep
	SCOREP_SWITCH=YES
    fi

    if [ ${PRODUCTION} -eq 1 ]
    then
	EXE_=${EXE}_Prod
	SCOREP_SWITCH=NO
    fi

    if [ ${ENERGY} -eq 1 ]
    then
	ENERGY_CPU_SWITCH=YES
	if [ "${GPU_SWITCH}" = "YES" ]
	then
	    ENERGY_GPU_SWITCH=YES
	else
	    ENERGY_GPU_SWITCH=NO
	fi
	EXE_=${EXE}_Energy
	SCOREP_SWITCH=NO
    else
	ENERGY_CPU_SWITCH=NO
	ENERGY_GPU_SWITCH=NO
    fi
    
    printf "\n\t Compiling ${EXE_SCOREP}... \n"
    COMPILE=compile_${EXE_}.txt
    make clean && make EXEC=${EXE_} COMPILER_CC=${CC} COMPILER_CPP=${CPP} DEBUG=${DEBUG_SWITCH} OMP=${OMP_SWITCH} FULL_GPU=${GPU_SWITCH} SYSTYPE=${SYSTEM} ENERGY_CPU=${ENERGY_CPU_SWITCH} ENERGY_GPU=${ENERGY_GPU_SWITCH} SCOREP=${SCOREP_SWITCH} |& tee ${COMPILE}
    file ${EXE_} |& tee -a ${COMPILE}
    ldd ${EXE_} |& tee -a ${COMPILE}

    for dir in ${OUT_DIR[@]}
    do
	cp ${EXE_} ${COMPILE} ${dir}
    done
    make clean
    rm -f ${EXE_} ${COMPILE}
done # loop over EXE

cd ${WORKDIR}

printf "\n\t Compilation done \n"
