#!/bin/bash
#SBATCH --nodes=                      
#SBATCH --ntasks-per-node=
#SBATCH --ntasks=
#SBATCH --cpus-per-task=
#SBATCH -p 
#SBATCH --time=00:00:00		       

###################################################################
#Script Name	: run-SN.sh    
#Description	: execute Dynemol from pwd using local dynemol src files 
#Args           : see usage       
#Author       	: Luis G C Rego
#date           : 04/Nov/2023
###################################################################

export DYNEMOLWORKDIR=$(pwd)
export DYNEMOLDIR=< path to dynemol directory >

##################################################
# environment commands here
##################################################

# Call antechamber.sh and pass all command-line arguments
$DYNEMOLDIR/antechamber.sh "$@"
source $DYNEMOLDIR/antechamber.sh

arguments=$($DYNEMOLDIR/antechamber.sh "$@")

mpirun $DYNEMOLDIR/dynemol $arguments
