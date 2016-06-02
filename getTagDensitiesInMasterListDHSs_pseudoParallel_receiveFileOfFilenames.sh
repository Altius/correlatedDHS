#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash getTagDensitiesInMasterListDHSs_pseudoParallel_receiveFileOfFilenames.sh min_idx max_idx

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ]; then
    echo -e "Usage:  $0 fileOfDensityFilenames minIdxInclusive maxIdxInclusive"
    echo -e "\twhere fileOfDensityFilenames contains the full path to each density file, one per line."
    exit 1
fi

fileOfFilenames=$1
if [ ! -s "$fileOfFilenames" ]; then
    echo -e "Error:  File $fileOfFilenames is empty, or was not found."
    exit 1
fi

ENCODE_AWG_BLACKLIST=wgEncodeHg19ConsensusSignalArtifactRegions.bed

MIN_IDX=$1
MAX_IDX=$2

OUTDIR=tagDensities
mkdir -p $OUTDIR

EXE=gchr
MASTER_DHSs=masterDHSs.bed3

#files=(density/A549-DS14289.75.20.uniques-density.36.hg19.bed.jarch
#    density/AG10803-DS12374.75.20.uniques-density.36.hg19.bed.jarch
#    density/AoAF-DS13513.75.20.uniques-density.36.hg19.bed.jarch
filenum=0
linenum=1
cat "$fileOfFilenames" | while read line
do
    if [ "$filenum" -lt "$MIN_IDX" ]; then
    ((filenum++))
    ((linenum++))
    continue
    fi
    if [ "$filenum" -gt "$MAX_IDX" ]; then
    break
    fi
    file=$line
    if [ ! -s "$file" ]; then
        echo "Error:  File $file,"
        echo "on line $linenum of $fileOfFilenames, was not found."
    exit 1
    fi
    is_jarch=$(basename "$file" | grep -c '\.jarch$')
    is_starch=$(basename "$file" | grep -c '\.starch$')
    if [ "$is_jarch" == "0" ] &&  [ "$is_starch" == "0" ]; then
    echo -e "Error:  Expect each density file to be in .jarch or .starch format. File $file,"
    echo -e "on line $linenum of $fileOfFilenames, does not end in .jarch or .starch."
    exit 1
    fi
    fnameStub=$(basename "$file" | cut -f1 -d '.')
    outfile=$OUTDIR/${fnameStub}_tagDensityInMasterDHSs.bed4
    # There are often many DHSs that are adjacent to each other,
    # such that one 20-bp density measurement can overlap 2 DHSs.
    # When 2 DHSs are overlapped, those will be printed in sort-bed order,
    # in a semicolon-delimited list.  Knowing they're in sorted order,
    # we can process them left-to-right in a fixed manner (i.e. without further comparison tests).
    # (Note: Each of those 2 DHSs gets a portion of the overlapping density element;
    # that portion is equal to the hand-calculated number of overlapping bp divided by 20bp;
    # we multiply by 0.05 instead of dividing by 20.)
    #
    # Each density value is the value for 150bp centered on the given 20bp.
    # Because density coordinates tile the genome beginning at chr<whatever>:65-85,
    # each DHS will be overlapped by exactly 7.5 20-bp regions (either 7 + .25 + .25 or 6 + .75 + .75).
    # So we divide the sum by 7.5 to get the mean density per 150bp.
    # (For efficiency's sake, we divide each density value by 7.5 before summing.)

    # It's exceedingly rare, and possibly nonexistent,
    # for a 20-bp density file entry to overlap a DHS boundary at one end
    # and the blacklist at the other end, so the assumption of 7.5 density entries per DHS
    # remains safe.  In the event that it does happen, something like 0.25/7.5 of the sum
    # will be zero.
    if [ "$is_jarch" == "1" ]; then
    "$EXE" "$file" \
        | bedops -n -1 - "$ENCODE_AWG_BLACKLIST" \
        | bedmap --count --echo --bases-uniq-f --echo-map - "$MASTER_DHSs" \
        | awk -F "|" '{split($2,x,"\t");dens=x[5]/7.5;
                         if(1==$1) {
                            print $4"\tid\t"dens*$3;
                         }
                         else if (2==$1) {
                            split($4,y,";");
                            split(y[1],z,"\t");
                            print y[1]"\tid\t"dens*(z[3]-x[2])*0.05;
                            split(y[2],z,"\t");
                            print y[2]"\tid\t"dens*(x[3]-z[2])*0.05;
                         }
                       }' \
            | bedmap --delim "\t" --prec 3 --echo --sum $MASTER_DHSs - \
        > "$outfile"
    else
    bedops -n -1 "$file" "$ENCODE_AWG_BLACKLIST" \
        | bedmap --count --echo --bases-uniq-f --echo-map - "$MASTER_DHSs" \
        | awk -F "|" '{split($2,x,"\t");dens=x[5]/7.5;
                         if(1==$1) {
                            print $4"\tid\t"dens*$3;
                         }
                         else if (2==$1) {
                            split($4,y,";");
                            split(y[1],z,"\t");
                            print y[1]"\tid\t"dens*(z[3]-x[2])*0.05;
                            split(y[2],z,"\t");
                            print y[2]"\tid\t"dens*(x[3]-z[2])*0.05;
                         }
                       }' \
            | bedmap --delim "\t" --prec 3 --echo --sum "$MASTER_DHSs" - \
        > "$outfile"
    fi
    ((filenum++))
    ((linenum++))
done

exit
