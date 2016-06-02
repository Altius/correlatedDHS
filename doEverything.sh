#! /bin/bash
# qsub -cwd -N JOB_ID -S doEverything.sh

set -e -o pipefail # quit when any section of any pipe, or single command, fails

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  $0 fileOfPeakfileNames.txt fileOfDensityfileNames.txt"
    exit 1
fi

peakfileNames=$1
densityfileNames=$2

qsub -cwd -N RunML -S run_master_list_simple__receiveInputFileOfFilenames.sh $peakfileNames
# Hmmm, need to check exit status....
