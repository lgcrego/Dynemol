# Stop on error
set -e

START_FOLDER=$(pwd)
RESULTS_FOLDER=$(pwd)/results

divider() {
    echo '=========='
}

copy_results() {
    # Copies the pre-defined files to the results folder
    #
    # Args:
    #   1. Destination folder path
    divider
    cp -v src/fort.13           $1
    cp -v src/frames-MM.pdb     $1
    cp -v src/MM_log.out        $1
    cp -v src/velocity_MM.out   $1
    divider
}

run_program() {
    # Executes the program saving all output to a file
    #
    # Args:
    #   1. File path to save program's output
    divider
    date
    divider
    time -p sh scripts/run.sh $2 $3 $4 2>&1 | tee $1
    divider
    date
    divider
}

battery_run() {
    # Create results folder for this test
    TEST_NAME=$1
    echo ''
    echo '=====>   START'
    echo "=====>   $TEST_NAME"
    echo ''

    mkdir -vp $RESULTS_FOLDER/$TEST_NAME

    # Enter program folder
    cd $TEST_NAME
    pwd

    # Sockets
    for s in 1 2
    do
        # Cores
        for c in 1 2 4 6 8 10
        do
            # Threads
            for t in 1 2
            do
                divider
                echo "Test s$s c$c t$t"
                divider

                # Create output folder
                TEST_NAME=$1
                TEST_SECTION="s$s""_c$c""_t$t"
                TEST_SECTION_FOLDER=$RESULTS_FOLDER/$TEST_NAME/$TEST_SECTION
                mkdir -vp $TEST_SECTION_FOLDER

                # Create output file
                OUT_FILE=$TEST_SECTION_FOLDER/out.txt
                run_program $OUT_FILE $s $c $t

                copy_results $TEST_SECTION_FOLDER

                divider
            done
        done
    done

    cd ..
    pwd

    echo ''
    echo "=====>   $TEST_NAME"
    echo '=====>   END'
    echo ''
}

# Create base results folder
mkdir -vp $RESULTS_FOLDER

battery_run $1
