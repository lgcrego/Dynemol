
# Compares the outputs of all the tests.
#
# It should not print anything if it is OK.
# Will print names of the files that are not equal.
find -name 'res*' -type d -exec \
        diff -rq -x velocity*.out -x out.txt res0 {} \;
