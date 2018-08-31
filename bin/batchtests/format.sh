# File that contains some scripts to be used to process the outputs from the tests.


format_text() {
    # Gets the text from a file and formats it to be processable by a sheet.
    #
    # Args:
    #   1: File to get text from.
    #
    # Output:
    # The final file format will be (without header and spaces):
    # Sockets   Cores   Threads     real    user    sys
    # 1         1       1           1.1     1.1     1.1
    # 1         1       1           1.1     1.1     1.1
    # ...

    # Gets the number of sockets, cores and threads from the folder names
    # and makes each number be in a new line
    REMOVED_FILE_NAMES=`perl -p -e 's/==>.*s(\d+)_c(\d+)_t(\d+).*/$1\n$2\n$3/' $1`

    # Removes 'real', 'user' and 'ssys' words
    REMOVED_WORDS=`echo "$REMOVED_FILE_NAMES" | perl -p -e 's/real|user|sys//'`

    # Removes all extra spaces on the file
    REMOVED_SPACES=`echo "$REMOVED_WORDS" | perl -p -e 's/ //'`

    # Replace dots for commas
    REMOVED_DOTS=`echo "$REMOVED_WORDS" | perl -p -e 's/\./,/'`

    # Make the file be 6 columns instead of a huge column
    echo "$REMOVED_DOTS" | paste - - - - - -
}


get_times() {
    # Gets the time of each out.txt file.
    #
    # The time is located in the last 3 lines of the file.
    #
    # Args:
    #   1: Which test to get the time of. Normally 'original'
    #       or 'protein'
    #
    # Out:
    #   {{ for every file }}
    #   ==> path/to/file.txt <==
    #   real ttt
    #   user ttt
    #   sys ttt
    find -name '*out.txt' -and -path "*$1*" -exec \
        tail -vn 3 {} \;
}


filter_by_socket() {
    # Filters the passed file by number of sockets.
    #
    # Args:
    #   1: Filename to filter
    #   2: Number of sockets
    #
    # Output:
    #   All lines from file that match the number of sockets.
    grep -E ^$2.* $1
}


get_times $1 > times_$1.txt

format_text times_$1.txt > format_$1.txt
rm times_$1.txt

filter_by_socket format_$1.txt $2
rm format_$1.txt