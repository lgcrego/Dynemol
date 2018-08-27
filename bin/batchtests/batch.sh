
# Runs the run.sh script the number of times passed as argument
# and creates a new folder for each run

set -e

for i in $(seq 0 `expr $1 - 1`)
do
    sh scripts/batchtests/test.sh original
    sh scripts/batchtests/test.sh protein
    mv -v results res$i
done


