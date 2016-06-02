#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash newCalcCorrelations_receiveThreshold.sh

source ~/.bashrc

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ] || [ "$5" == "" ] || [ "$6" == "" ]; then
    echo -e "Usage:  $0 infile THRESHOLD outfile_all_promotersFirst.bed8 outfile_aboveTHRESHOLD_promotersFirst.bed8 outfile_all_distalsFirst.bed8 outfile_aboveTHRESHOLD_distalsFirst.bed8"
    echo -e "\twhere THRESHOLD is traditionally 0.7, currently 0.5"
    exit
fi

INFILE=$1
THRESHOLD=$2
OUTFILE_ALL_PROM_FIRST=$3
OUTFILE_THRESHOLDED_PROM_FIRST=$4
OUTFILE_ALL_DISTAL_FIRST=$5
OUTFILE_THRESHOLDED_DISTAL_FIRST=$6

EXE=/net/lebowski/vol1/work/erynes/data/DNasePatternMatching/gencode10/24Jul13andBeyond/33celltypes/CalcCorrelations_30Jul13

$EXE $INFILE $OUTFILE_ALL_PROM_FIRST

if [ -s "$OUTFILE_ALL_PROM_FIRST" ]; then
    awk -v threshold=$THRESHOLD '{if($8>=threshold){print;}}' $OUTFILE_ALL_PROM_FIRST > $OUTFILE_THRESHOLDED_PROM_FIRST
    awk 'BEGIN{OFS="\t"}{print $5,$6,$7,$4,$1,$2,$3,$8;}' $OUTFILE_ALL_PROM_FIRST \
	| sort-bed - \
	> $OUTFILE_ALL_DISTAL_FIRST
    awk -v threshold=$THRESHOLD '{if($8>=threshold){print;}}' $OUTFILE_ALL_DISTAL_FIRST > $OUTFILE_THRESHOLDED_DISTAL_FIRST
fi

exit
