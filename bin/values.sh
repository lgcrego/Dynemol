# File created to get values from some file based on the text that comes before it.
#
# For example, if I have a file, test.txt, with the contents:
#
#   time ==> 0
#   time ==> 1
#   asdasdasdas
#
# I can get the values 0 and 1 with the command: "sh values.sh test.txt 'time ==>'"

# Pattern to check for
PATTERN=$2

# All lines that have the pattern
LINES=$(grep "$PATTERN" $1)

# Remove pattern from the lines and output the values only
echo "$LINES" | sed "s/$PATTERN//g"