
# Compare the files from src with the files from benchmark folder

help() {
    echo ''
    echo 'Compare current outputs of the src folder with desired benchmark'
    echo 'folder'
    echo ''
    echo 'Args:'
    echo '  1. Folder, in benchmark/ to compare files from'
    echo ''
    echo 'Output:'
    echo '  Same as diff -q output for every file'
    echo ''
}

# Check for args
if [ $# -ne 1 ]
then
    help
    exit
fi

# Compare
diff -q ./benchmark/$1/fort.13          src/fort.13
diff -q ./benchmark/$1/frames-MM.pdb    src/frames-MM.pdb
diff -q ./benchmark/$1/MM_log.out       src/MM_log.out
diff -q ./benchmark/$1/velocity_MM.out  src/velocity_MM.out
