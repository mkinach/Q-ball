#! /bin/bash

#=======================================================================
# Functions
#=======================================================================

# generate a complete PAMR parameter file
genparam() {
  cat "${INPUT_PATH}/common.fparam" "${INPUT_PATH}/${1}.fparam" \
      "${INPUT_PATH}/initdata.fparam" "${INPUT_PATH}/common.rtparam" \
      "${INPUT_PATH}/${1}.rtparam" > "${OUTPUT_PATH}/pamr.allparam"
}

# print usage instructions and exit
usage() {
  echo
  echo "usage: ${SCRIPT_NAME} [-m <mode>] [-n <integer>] [-v <potential>]"
  echo
  echo "    -m <mode>       Specify a PAMR mode to use ('mg' or 'evo' supported; default is 'evo')"
  echo "    -n <integer>    Specify the number of cores to use (integer > 0 supported; default is 1)"
  echo "    -v <potential>  Generate initial data corresponding to a specific scalar potential ('log' or 'poly' supported)"
  echo "    -h              Print this help message"
  echo
  exit 1
}

#=======================================================================
# Main Program
#=======================================================================

# get the script name
SCRIPT_NAME=$(basename "$0")

# define relevant path variables
PAMR_PATH=$(realpath ./src/pamr)
INPUT_PATH=$(realpath ./run/input)
OUTPUT_PATH=$(realpath ./run/output)

# default settings to use
NUMCORES=1  # number of cores
MODE=evo    # which PAMR "mode" to use (multigrid or time evolution)

# parse command line options, if any
while getopts "n:m:v:h" OPT; do
case ${OPT} in
  n)
    NUMCORES=${OPTARG}
    if ! [[ ${NUMCORES} =~ ^[1-9][0-9]*$ ]]; then usage; fi
    ;;
  m)
    MODE=${OPTARG}
    if [[ ${MODE} != "mg" && ${MODE} != "evo" ]]; then usage; fi
    ;;
  v)
    POTENTIAL=${OPTARG}
    if [[ ${POTENTIAL} != "log" && ${POTENTIAL} != "poly" ]]; then usage; fi
    ;;
  h) usage ;;
  \?) usage ;;
  esac
done

# check for invalid options
if [[ ${MODE} == "mg" && ${NUMCORES} -gt "1" ]]; then
  printf "\nMODE=mg not supported with NUMCORES>1\n\n" && exit 1
fi

# clean up the directories
make vclean

# generate initial data
if [ -n "$POTENTIAL" ]; then
  (
    cd "${INPUT_PATH}" || exit 1
    maple "shoot_${POTENTIAL}.mpl" \
      && printf "\nInitial data generation complete." \
      && printf " Please verify before proceeding.\n"
  )
elif [ ! -e "${INPUT_PATH}/initdata" ]; then
  printf "\nInitial data file not found.\n"
  usage
fi

# generate parameter file for evolution
genparam evo
printf "\nPreparing to run with MODE=%s and NUMCORES=%s\n\n" "${MODE}" "${NUMCORES}"
genparam "${MODE}"  # regenerate the multigrid parameter file (if necessary)

# check whether PAMR debug options are enabled
DEBUG=$(cat "${PAMR_PATH}"/Makefile | grep -c '^FDEBUGFLAGS')

# confirm settings before compilation
read -r -p 'Proceed? (Y or n): ' CHOICE
case "${CHOICE}" in
  Y) ;;
  *) exit 0;;
esac

# compile RNPL and PAMR code
printf "\n\n"
make all
printf "\n\n"

# run the program
cd "${OUTPUT_PATH}" || exit 1
if [[ ${DEBUG} -gt 0 ]]; then
  printf "\n\e[01;41m DEBUG FLAG ON! \e[0m \n\n"
  printf "You should enter 'run pamr.allparam' into the gdb instance(s)\n\n"
  mpirun --oversubscribe -np "${NUMCORES}" konsole -e gdb "${PAMR_PATH}/qball"
else
  mpirun --oversubscribe -np "${NUMCORES}" "${PAMR_PATH}/qball" pamr.allparam
fi

# if running in multigrid mode, automatically proceed with the time evolution
if [[ ${MODE} == "mg" ]]; then
  MODE=evo
  genparam ${MODE}
  printf "Running with MODE=%s, NUMCORES=%s\n\n" "${MODE}" "${NUMCORES}"
  if [[ ${DEBUG} -gt 0 ]]; then
    printf "\n\e[01;41m DEBUG FLAG ON! \e[0m \n\n"
    printf "You should enter 'run pamr.allparam' into the gdb instance(s)\n\n"
    mpirun --oversubscribe -np "${NUMCORES}" konsole -e gdb "${PAMR_PATH}/qball"
  else
    mpirun --oversubscribe -np "${NUMCORES}" "${PAMR_PATH}/qball" pamr.allparam
  fi
fi
