#!/bin/bash

# get the executable
cd ${OUT_DIR}
EXEC=($(find $(realpath ./) -maxdepth 1 -name "*Scorep" -executable -type f -print))
if [ ${#EXEC[@]} -lt 1 ]
then
    printf "\n\t Cannot find any Scorep executable... aborting... \n"
    exit 1
fi

export SCOREP_VERBOSE=1
export SCOREP_ENABLE_PROFILING=0
export SCOREP_ENABLE_TRACING=1
export SCOREP_EXPERIMENT_DIRECTORY=tracing

export OMP_WAIT_POLICY=ACTIVE

for EXE in ${EXEC[@]}
do
    printf "\n\t Running ${EXE}... \n"

    PROFILING=($(find $(realpath ./) -type d -name "$(basename ${EXE})_profiling" -print))
    if [ ${#PROFILING[@]} -eq 0 ]
    then
	printf "\n\t 'profiling' folder not found for ${EXE}... aborting... \n\n"
	exit 2
    elif [ ${#PROFILING[@]} -ge 2 ]
    then
	printf "\n\t more than one 'profiling' folder found for ${EXE}... aborting... \n\n"
	exit 4
    else
	cd ${PROFILING[0]}

	for NT in ${NTASKS[@]}
	do
	    NTASKS_PER_NODE=$((NT / NODES))

	    for MAP in ${MPI_MAP_BY[@]}
	    do
		case ${MAP} in

		    socket)
			PPR=$((NT / SOCKETS))
			if [[ ${PPR} -eq 0 ]]
			then
			    PPR=1
			fi
			;;

		    numa)
			PPR=$((NT / NUMA))
			if [[ ${PPR} -eq 0 ]]
			then
			    PPR=1
			fi
			;;
		    *)
			printf "\n\t Unknown \n"
			exit 7
			;;	
		esac

    		for OMP in ${OMP_THR[@]}
    		do
    		    # set the number of omp threads
    		    export OMP_NUM_THREADS=${OMP}

		    # get the subdirectory each MAP/MPI/OMP configuration
		    SUBDIR=nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}
		    SUB_DIR=($(find $(realpath ./) -maxdepth 1 -name "${SUBDIR}" -type d -print))
		    if [ ${#SUB_DIR[@]} -eq 0 ]
		    then
			printf "\n\t Cannot find ${PWD}/${SUBDIR}... aborting... \n\n"
			exit 4
		    fi
		    
		    # enter the directory
		    cd ${SUB_DIR[0]}
		    
		    # get the profiling folder
		    PFOLDER=($(find $(realpath ./) -maxdepth 1 -name "profiling" -type d -print))
		    if [ ${#PFOLDER[@]} -eq 0 ]
		    then
			printf "\n\t Cannot find ${PWD}/profiling... aborting... \n\n"
			exit 5
		    fi

		    FILTERING=$(find $(realpath ./profiling) -maxdepth 1 -name "custom_scorep.filter" -type f -print)
		    if [ ! -f ${FILTERING} ]
		    then
			printf "\n\t Cannot find 'custom_scorep.filter' ... aborting... \n\n"
			exit 6
		    else
			export SCOREP_FILTERING_FILE=${FILTERING}
		    fi
		    
		    MOVE_TO=${PFOLDER[0]}_tracing
		    OUT=${PWD}/$(basename ${EXE})_nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}_tracing_output.txt
		    
		    printf "\n\t Running ${EXE} using:"                                            |& tee ${OUT}
		    printf "\n\t                       ${NODES} nodes"                             |& tee -a ${OUT}
		    printf "\n\t                       ${NT} MPI processes"                        |& tee -a ${OUT}
		    printf "\n\t                       ${NTASKS_PER_NODE} MPI processes per node"  |& tee -a ${OUT}
		    printf "\n\t                       MPI processes mapped by ${MAP}"             |& tee -a ${OUT}
		    printf "\n\t                       each MPI process spawns ${OMP} OMP threads" |& tee -a ${OUT}
		    printf "\n\t                       $((NT * OMP)) processors used\n\n"          |& tee -a ${OUT}
		    printf "\n\t mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARAMFILE} \n" |& tee -a ${OUT}
		    
		    mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARAMFILE} |& tee -a ${OUT}

		    # change directory name
		    mv ${PFOLDER[0]} ${MOVE_TO}
		    mv ${SCOREP_EXPERIMENT_DIRECTORY}/* ${MOVE_TO}
		    rm -rf ${SCOREP_EXPERIMENT_DIRECTORY}

		    cd ../
    		done # OMP in OMP_THREADS
	    done # MAP in MAP_BY
	done # NT in NTASKS

	cd ..
	
    fi # profiling folder
done

cd ${WORKDIR}

printf "\n\t Tracing done \n"
