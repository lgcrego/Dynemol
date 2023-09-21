#!/bin/bash 

###################################################################
#Script Name	: run-SN.sh                                                                                             
#Description	: execute Dynemol from pwd using local dynemol src files                                                                               
#Args           : see usage                                                                                          
#Author       	: Luis G C Rego
#date           : 21/Sept/2023
###################################################################
#!/bin/bash


# Function to display usage instructions
usage() {
  echo "Usage: $0 [OPTIONS]"
  echo " "
  echo "  -h, --help, help     Display this help message"
  echo " "
  echo "  resume,   resumes MD simulation; need to mv velocities.out -> velocities.inpt; no need to update input.pdb"
  echo " "
  echo "  preview,  just initiates Dynemol, for the sake of checking initial execution"
  echo " "
  echo "  spawn     generates thermal configurations out of an MD run ; usage below"
  echo "  spawn     {# of configurations to be saved}"
  echo " "
  echo "  PDOS      generates PDOS of atomic features , implemented in driver = diagnostic ; usage below"
  echo "  PDOS      {EHSymbol or Symbol}"
  echo " "
  echo "  MO        generates MO cube files for visualization , implemented in driver = {diagnostic,Genetic_Alg} ; usage below"
  echo "  MO        {MO_i}-{MO_f}             ,  range of MO's from MO_1 to MO_f"
  echo "  MO        {MO_1} {MO_2} ... {MO_n}  ,  list of MO's separated by spaces"
  echo " "
  echo "  for (driver = MM_Dynamics) and (driver_MM = Parametrize) choose among the arguments below"
  echo "  newOPT , repeat , resume"
  exit 1
}

# Parse command line arguments
case "$1" in
  resume)
    arguments="resume" 
    ;;
  preview)
    arguments="preview" 
    ;;
  spawn)
    arguments="spawn $2" 
    ;;
  PDOS|pdos|Pdos)
    if [ "$#" -lt 2 ]; then
         usage; exit
    fi
    arguments="PDOS $2" 
    ;;
  -h|--help|help)
    usage
    ;;
  *)
    echo "Error: Unknown option: $1"
    usage
    ;;
esac



export DYNEMOLWORKDIR=$(pwd)
export DYNEMOLDIR=<path to dynemol dir>

$DYNEMOLDIR/dynemol $arguments

