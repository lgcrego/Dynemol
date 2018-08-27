#!/bin/bash

# Stops when error is found
set -e

# Prints help message
help() {
    echo ''
    echo 'Runs the project with the desired, thread affinity related, enviroment'
    echo 'variables set'
    echo ''
    echo '  - If no 4th argument has been provided, it will be set to 1'
    echo ''
    echo 'Args:'
    echo '  1. Number of sockets'
    echo '  2. Number of cores per socket'
    echo '  3. Number of threads per core'
    echo '  4. Number of MKL threads'
    echo ''
}


run() {
    # Enviroment vaiables used to test
    export KMP_HW_SUBSET=${1}s,${2}c,${3}t
    export MKL_NUM_THREADS=$4
    export KMP_AFFINITY='verbose,granularity=fine,compact'


    # Run (it must be ran from inside the folder)
    cd src
    ./a
}


# 3 args
if [ $# -eq 3 ]
then
    run $1 $2 $3 1
    exit
fi

if [ $# -eq 4 ]
then
    run $1 $2 $3 $4
    exit
fi

help
