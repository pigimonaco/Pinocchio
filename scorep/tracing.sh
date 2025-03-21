#!/bin/bash

for ((index=0 ; index<${#OUT_DIR[@]} ; index++))
do
    DIRECTORY=${OUT_DIR[index]} ; cd ${DIRECTORY}

    EXEC=($(find $(realpath ./) -maxdepth 1 -name "*_Scorep" -executable -type f -print))
    if [ ${#EXEC[@]} -lt 1 ]
    then
	printf "\n\t Cannot find any _Scorep executable... aborting... \n"
	exit 1
    fi

    export SCOREP_VERBOSE=1
    export SCOREP_ENABLE_PROFILING=0
    export SCOREP_ENABLE_TRACING=1
    export SCOREP_METRIC_PAPI= # empty
    export SCOREP_METRIC_PERF= # empty
    export SCOREP_EXPERIMENT_DIRECTORY=tracing

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
	    exit 3
	else
	    # enter the profiling directory
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
			    exit 4
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
			    exit 5
			fi
			
			# enter the directory
			cd ${SUB_DIR[0]}
			
			# get the profiling folder(s) (more than one if PAPI events were collected)
			PFOLDER=($(find $(realpath ./) -maxdepth 2 -name "profiling" -type d -print))
			if [ ${#PFOLDER[@]} -eq 0 ]
			then
			    printf "\n\t Cannot find ${PWD}/profiling... aborting... \n\n"
			    exit 6
			fi

			TARGET="/profiling"
			# loop over profiling folders
			for PROF in ${PFOLDER[@]}
			do
			    cd ${PROF::-${#TARGET}}
			    FILTERING=$(find $(realpath ./profiling) -maxdepth 1 -name "custom_scorep.filter" -type f -print)
			    if [ ! -f ${FILTERING} ]
			    then
				printf "\n\t Cannot find 'custom_scorep.filter' ... aborting... \n\n"
				exit 7
			    else
				export SCOREP_FILTERING_FILE=${FILTERING}
			    fi
			
			    PARAMFILE=($(find ./ -maxdepth 1 -name "paramfile_*" -type f -print))
			    if [ ${#PARAMFILE[@]} -ne 1 ]
			    then
				printf "\n\t Only one paramfile is expected in ${PWD} ...aborting... \n"
				exit 8
			    fi
			
			    OUT=${PWD}/$(basename ${EXE})_nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}_tracing_output.txt
			
			    STREAM="\n\t Running ${EXE} using:"
			    STREAM+="\n\t                       ${NODES} nodes"
			    STREAM+="\n\t                       ${NT} MPI processes"
			    STREAM+="\n\t                       ${NTASKS_PER_NODE} MPI processes per node"
			    STREAM+="\n\t                       MPI processes mapped by ${MAP}"
			    STREAM+="\n\t                       each MPI process spawns ${OMP} OMP threads"
			    STREAM+="\n\t                       $((NT * OMP)) processors used\n\n"
			    STREAM+="\n\t mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARAMFILE[0]} \n\n"
			    echo -e ${STREAM} |& tee ${OUT}

			    mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARAMFILE[0]} |& tee -a ${OUT}

			    # change directory name
			    MOVE_TO=${PROF}_tracing
			    mv ${PROF} ${MOVE_TO}
			    mv ${SCOREP_EXPERIMENT_DIRECTORY}/* ${MOVE_TO}
			    rm -rf ${SCOREP_EXPERIMENT_DIRECTORY}

			    cd -
			done # loop over PFOLDER
			cd ../
    		    done # OMP in OMP_THREADS
		done # MAP in MAP_BY
	    done # NT in NTASKS
	    cd ${DIRECTORY}
	fi # profiling folder
    done # loop over EXE
done # loop over OUT_DIR

cd ${WORKDIR}

printf "\n\t Tracing done \n"
