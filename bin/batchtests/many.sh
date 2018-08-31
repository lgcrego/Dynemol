
# Runs the run.sh script the number of times passed as argument
# and creates a new folder for each run

set -e

for i in $(seq 0 `expr $1 - 1`)
do
    find -name run.sh -exec sh {} protein original \;
    mv -v results res$i
done


