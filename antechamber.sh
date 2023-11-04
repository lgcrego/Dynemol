#!/bin/bash 

###################################################################
#Script Name	: run-SN.sh                                                                                             
#Description	: execute Dynemol from pwd using local dynemol src files                                                                               
#Args           : see usage                                                                                          
#Author       	: Luis G C Rego
#date           : 30/Sept/2023
###################################################################

# Function to display usage instructions
usage() {
  echo "Usage: $0 [OPTIONS]"
  echo " "
  echo "  -h, --help, help     Display this help message"
  echo " "
  echo "  manipulate,   calls suite of manipulate routines for pre- and post-processing data"
  echo " "
  echo "  resume,       resumes MD simulation; need to mv velocities.out -> velocities.inpt; no need to update input.pdb"
  echo " "
  echo "  preview,      just initiates Dynemol, for the sake of checking initial execution"
  echo " "
  echo "  spawn,        generates thermal configurations out of an MD run ; usage below"
  echo "  usage:        spawn {# of configurations to be saved}"
  echo " "
  echo "  PDOS,         generates PDOS of atomic features , implemented in driver = diagnostic ; usage below"
  echo "  usage:"
  echo "  PDOS {EHSymbol or Symbol}"
  echo "  PDOS {AO} {MO_i} - {MO_f} ,  range of MO's from MO_i to MO_f"
  echo " "
  echo "  MO,           generates MO cube files for visualization , implemented in driver = {diagnostic,Genetic_Alg} ; usage below"
  echo "  usage:"
  echo "  MO   {MO_i} - {MO_f}      ,  range of MO's from MO_i to MO_f"
  echo " "
  echo "  for (driver = MM_Dynamics) and (driver_MM = Parametrize) choose among the arguments below"
  echo "  newOPT , repeat , resume"
  exit 1
}

force_flag=false

# Parse command line arguments
if [[ "$#" -ne 0 ]]; then
    case "$1" in
      -h|--help|help)
        usage
        ;;
      -f|--force)
        force_flag=true
        shift
        ;;
      manipulate)
        $DYNEMOLDIR/manipulate/manipulate
        exit 1
        ;;
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
        elif [ "$#" -eq 2 ]; then 
             arguments="PDOS $2" 
        elif [ "$#" -eq 5 ]; then 
             arguments="PDOS AO $3"\ -\ "$5" 
        fi
        ;;
      MO|mo)
        if [ "$#" -le 2 ]; then
             usage; exit
        elif [ "$#" -eq 3 ]; then 
             arguments="MO $2 $3" 
        elif [ "$#" -eq 4 ]; then 
             arguments="MO $2"\ -\ "$4" 
        fi
        ;;
      *)
        echo "Error: Unknown option: $1"
        usage
        ;;
    esac
fi
echo "$arguments"
export force_flag
