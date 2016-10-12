#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash makeMasterFileWithTagSumVectors.sh

PID=$$
TEMP1=tempfile1_${PID}.bed
TEMP2=tempfile2_${PID}.bed
outfile=masterDHSsAndTagCounts.bed4

TT=totalTagCountsPerSample.txt

files=(tagDensities/*)

MAX_COUNT=0
set -x
for file in ${files[*]}
do
    celltype=$(basename "$file" _tagDensityInMasterDHSs.bed4)
    count=$(grep -w "$celltype" $TT | cut -f2)
    if [ "$count" -gt "$MAX_COUNT" ]; then
        MAX_COUNT=$count
    fi
done
if [ "$MAX_COUNT" == "0" ]; then
    echo -e "Error:  Failed to find maximum genomewide tag count. Exiting...."
    exit
fi

file="${files[0]}"
celltype=$(basename "$file" _tagDensityInMasterDHSs.bed4)
count=$(grep -w "$celltype" "$TT" | cut -f2)
if [ "$count" == "" ]; then
    echo -e "Failed to find $celltype in $TT."
    exit
fi

awk -v count="$count" -v maxCount="$MAX_COUNT" 'BEGIN{OFS="\t"}{print $1,$2,$3,int(($4*1.0*maxCount/count) + 0.5);}' "$file" > "$TEMP2"

for file in ${files[*]}
do
    if [ "$file" == "${files[0]}" ]; then
        continue
    fi
    celltype=$(basename "$file" _tagDensityInMasterDHSs.bed4)
    count=$(grep -w "$celltype" "$TT" | cut -f2)
    if [ "$count" == "" ]; then
        echo -e "Failed to find $celltype in $TT."
        exit
    fi
    awk -v count="$count" -v maxCount="$MAX_COUNT" 'BEGIN{OFS="\t"}{print int(($4*1.0*maxCount/count) + 0.5);}' "$file" \
        | paste -d "," $TEMP2 - \
        > $TEMP1
    mv $TEMP1 $TEMP2
done
mv $TEMP2 $outfile
rm -f $TEMP1
exit
