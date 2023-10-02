#!/bin/bash 

###################################################################
#Script Name	: run-SN.sh                                                                                             
#Description	: execute Dynemol from pwd using local dynemol src files                                                                               
#Args           : see usage                                                                                          
#Author       	: Luis G C Rego
#date           : 30/Sept/2023
###################################################################

set -e  # Enable the "exit on error" option

export DYNEMOLWORKDIR=$(pwd)
export DYNEMOLDIR=< path to dynemol directory >

# Call antechamber.sh and pass all command-line arguments
$DYNEMOLDIR/antechamber.sh "$@"

arguments=$($DYNEMOLDIR/antechamber.sh "$@")

$DYNEMOLDIR/dynemol $arguments

