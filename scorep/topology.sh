#!/bin/bash

# get the hw resources using the 'lscpu' tool

if [ $# -lt 1 ]
then
    printf "\n\t Usage: ./topology.sh <number of computational nodes> \n\n"
    exit 1
else
    NNODES=$1
fi

TOPOLOGY=${HOSTNAME}_topology.txt
lscpu |& tee ${TOPOLOGY}
# CPU model
CPU=$(lscpu        | grep "Model name"   | gawk -F ":" '{print $2}' | xargs)
# hwthreads per core
HW_THREADS=$(lscpu | grep "Thread(s)"    | gawk -F ":" '{print $2}' | xargs)
# hw cores per socket 
HW_CORES=$(lscpu   | grep "Core(s)"      | gawk -F ":" '{print $2}' | xargs)
# sockets
HW_SOCKETS=$(lscpu | grep "Socket(s)"    | gawk -F ":" '{print $2}' | xargs)
# numa nodes
HW_NUMA=$(lscpu    | grep "NUMA node(s)" | gawk -F ":" '{print $2}' | xargs)
# L1d cache
CACHE_L1D=$(lscpu  | grep "L1d"          | gawk -F ":" '{print $2}' | xargs)
# L1i cache
CACHE_L1I=$(lscpu  | grep "L1i"          | gawk -F ":" '{print $2}' | xargs)
# L2 cache
CACHE_L2=$(lscpu   | grep "L2"           | gawk -F ":" '{print $2}' | xargs)
# L3 cache
CACHE_L3=$(lscpu   | grep "L3"           | gawk -F ":" '{print $2}' | xargs)

# get the overall hw available resources
THREADS=$((HW_THREADS * HW_CORES   * HW_SOCKETS * NNODES)) # total number of available hwthreads
CORES=$((  HW_CORES   * HW_SOCKETS              * NNODES)) # total number of available cores
SOCKETS=$((HW_SOCKETS                           * NNODES)) # total number of available sockets
NUMA=$((   HW_NUMA                              * NNODES)) # total number of available numa regions

printf "\n\t CPU architecture:"
printf "\n\t\t Model name         : ${CPU}"
printf "\n\t\t Thread(s) per core : ${HW_THREADS}"
printf "\n\t\t Core(s) per socket : ${HW_CORES}"
printf "\n\t\t Socket(s)          : ${HW_SOCKETS}"
printf "\n\t\t NUMA node(s)       : ${HW_NUMA}"
printf "\n\t Caches:"
printf "\n\t\t L1d                : ${CACHE_L1D}"
printf "\n\t\t L1i                : ${CACHE_L1I}"
printf "\n\t\t L2                 : ${CACHE_L2}"
printf "\n\t\t L3                 : ${CACHE_L3}"
printf "\n\n\t NODES: ${NNODES}"
printf "\n\t\t THREADS: ${THREADS} - CORES: ${CORES} - SOCKETS: ${SOCKETS} - NUMA: ${NUMA} \n\n"
