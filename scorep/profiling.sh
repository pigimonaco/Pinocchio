#!/bin/bash

function generate_filter_file
{
    # generating initial filter file
    if [ -f ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex ]
    then
	# performance analysis    
	ANALYSIS=${PWD}/$(basename ${EXE})_nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}_scorep_profiling_analysis.txt
	# create the initial_scorep.filter
	scorep-score -r -m -g -s totaltime ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex |& tee ${ANALYSIS}

	# write a custom filter including COM group only
	_FILTER=custom_scorep.filter
	touch ${_FILTER}
	printf "SCOREP_REGION_NAMES_BEGIN\n" >> ${_FILTER}
	printf "EXCLUDE *\n" >> ${_FILTER}
	COM_group=$(grep -w "COM" ${ANALYSIS} | gawk 'NR>1 {print $7}' | xargs)
	for item in ${COM_group}
	do
            printf "INCLUDE ${item}\n" >> ${_FILTER}
	done
	printf "SCOREP_REGION_NAMES_END\n" >> ${_FILTER}

	# test the effect of the initial filter
	FILTER=$(find $(realpath ./) -name "initial_scorep.filter" -type f -print)
	if [ -f ${FILTER} ]
	then
	    # analyze the effect of the initial filter         
	    printf "\n\t Analyzing the effect of the initial filter \n\n" |& tee -a ${ANALYSIS}
	    scorep-score -r -m -s totaltime -f ${FILTER} ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex |& tee -a ${ANALYSIS}
	    mv ${FILTER} ${SCOREP_EXPERIMENT_DIRECTORY}
	else
	    printf "\n\t Cannot generate score-p initial filter file \n" |& tee -a ${ANALYSIS}
	fi # FILTER

	# test the effect of the custom filter
	printf "\n\t Analyzing the effect of the custom filter \n\n" | tee -a ${ANALYSIS}
	scorep-score -r -m -s totaltime -f ${_FILTER} ${SCOREP_EXPERIMENT_DIRECTORY}/profile.cubex |& tee -a ${ANALYSIS}
	mv ${_FILTER} ${ANALYSIS} ${SCOREP_EXPERIMENT_DIRECTORY}
    fi
}

for ((index=0 ; index<${#OUT_DIR[@]} ; index++))
do
    DIRECTORY=${OUT_DIR[index]} ; cd ${DIRECTORY}

    EXEC=($(find . -name "*_Scorep" -maxdepth 1 -executable -type f -print))
    if [ ${#EXEC[@]} -lt 1 ]
    then
	printf "\n\t Cannot find any _Scorep executable... aborting... \n"
	exit 1
    fi

    PARAMETER_FILE=($(find . -maxdepth 1 -name "paramfile_*" -type f -print))
    if [ ${#PARAMETER_FILE[@]} -ne 1 ]
    then
	printf "\n\t Only one paramfile is expected in ${DIRECTORY} ...aborting... \n"
    fi

    PARAMFILE_=${PWD}/$(basename ${PARAMETER_FILE[0]})
    
    export SCOREP_VERBOSE=1
    export SCOREP_ENABLE_PROFILING=1
    export SCOREP_ENABLE_TRACING=0
    # export SCOREP_METRIC_PERF=${PERF_METRIC}
    export SCOREP_EXPERIMENT_DIRECTORY=profiling

    # store the output files
    FILE=

    # loop over executables
    for EXE in ${EXEC[@]}
    do
	printf "\n\t Running ${EXE}... \n"

	EXE=${PWD}/$(basename ${EXE})

	# create one directory for each executable
	DIR_NAME=${PWD}/$(basename ${EXE})_${SCOREP_EXPERIMENT_DIRECTORY}
	mkdir -p ${DIR_NAME} && cd ${DIR_NAME}

	# loop over MPI tasks
	for NT in ${NTASKS[@]}
	do
	    NTASKS_PER_NODE=$((NT / NODES))

	    # loop over MPI mapping
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
			exit 2
			;;	
		esac

		# loop over OMP threads
    		for OMP in ${OMP_THR[@]}
    		do
    		    # set the number of omp threads
    		    export OMP_NUM_THREADS=${OMP}

    		    # create a subdirectory for each MAP/MPI/OMP configuration
		    SUB_DIR=${PWD}/nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}
		    mkdir -p ${SUB_DIR} && cd ${SUB_DIR} && rm -rf * && cp ${OUTPUTS} .

		    PARA=$(basename ${PARAMFILE_})
		    # set the MaxMem parameter according to NTASKS_PER_NODE
		    MAX_MEM_MPI_TASK=$(( MAX_MEM * 1024 / NTASKS_PER_NODE )) # MaxMem in Mbyte
		    gawk -v max_mem_mpi_task=${MAX_MEM_MPI_TASK} '{if ($1=="MaxMem") { $2 = max_mem_mpi_task } print $0}' ${PARAMFILE_} > tmp.txt && mv tmp.txt ${PARA}

		    STREAM="\n\t Running ${EXE} using:"
		    STREAM+="\n\t                       ${NODES} nodes"
		    STREAM+="\n\t                       ${NT} MPI processes"
		    STREAM+="\n\t                       ${NTASKS_PER_NODE} MPI processes per node"
		    STREAM+="\n\t                       MPI processes mapped by ${MAP}"
		    STREAM+="\n\t                       each MPI process spawns ${OMP} OMP threads"
		    STREAM+="\n\t                       $((NT * OMP)) processors used"
		    STREAM+="\n\t mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARA}"

		    if [ ${#PAPI_METRIC} -ne 0 ]
		    then
			# loop over PAPIs
			for ((papi_index=0 ; papi_index<${#PAPI_METRICS[@]} ; papi_index++))
			do
			    mkdir -p ${PAPI_OUT_DIR[papi_index]} && cd ${PAPI_OUT_DIR[papi_index]} && cp ${OUTPUTS} ../${PARA} .
    			    OUT=${PWD}/$(basename ${EXE})_nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}_${PAPI_METRICS[papi_index]}_profiling_output.txt
    			    FILE+=("${OUT}")
			    export SCOREP_METRIC_PAPI=${PAPI_METRICS[papi_index]}
			    STREAM+="\n\t                       PAPI_METRIC: ${SCOREP_METRIC_PAPI}\n\n"
			    # print to file
			    echo -e ${STREAM} |& tee ${OUT}
			    # timer
    			    SECONDS=0
			    # run the code
			    mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARA} |& tee -a ${OUT}
    			    # get execution time
    			    TIME=$(($SECONDS))
    			    printf "\n\n\t Execution time is ${TIME} seconds \n\n" |& tee -a ${OUT}
			    # generating initial filter file
			    generate_filter_file
			    cd -
			done # loop papi_index
		    else
			OUT=${PWD}/$(basename ${EXE})_nodes_${NODES}_map_${MAP}_MPI_${NT}_OMP_${OMP}_profiling_output.txt
    			FILE+=("${OUT}")
			# print to file
			echo -e ${STREAM} |& tee ${OUT}
			# timer
    			SECONDS=0
			export SCOREP_METRIC_PAPI=
			mpirun -n ${NT} --map-by ppr:${PPR}:${MAP}:PE=${OMP} --bind-to core --report-bindings ${EXE} ${PARA} |& tee -a ${OUT}
    			# get execution time
    			TIME=$(($SECONDS))
    			printf "\n\n\t Execution time is ${TIME} seconds \n\n" |& tee -a ${OUT}
			# generating initial filter file
			generate_filter_file
		    fi

		    cd ../
    		done # OMP in OMP_THREADS
	    done # MAP in MAP_BY
	done # NT in NTASKS

	cd ${DIRECTORY}
    done # EXE in EXEC

    for TXT in ${FILE[@]}
    do
	RET=$(grep -i "Execution " ${TXT})
	printf "\n\t ${TXT}: ${RET}"
    done
done # loop over OUT_DIR

cd ${WORKDIR}

printf "\n\t Profiling done \n"
