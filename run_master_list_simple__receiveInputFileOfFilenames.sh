#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash run_master_list_simple__receiveInputFileOfFilenames.sh

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  $0 sampleFile.txt outfilename [blasklistfile]"
    echo -e "\twhere sampleFile contains the full path to each peak file in column 1, one per line."
    echo -e "\tblacklist file is optional"
    exit 1
fi

fileOfFilenames=$1

if [ "$3" == "" ]; then
    BLACKLIST=/dev/null
else
    BLACKLIST=$3
fi

## Make master set of unique DHS positions in the genome.  This will be the union of
## the MCV peaks and all single-cell strict peaks not in MCV zones.
## Also extract the "cell-type predominant" set.  This is everything that's left over from the
## single-cell line union after you subtract off MCV and CTS.
out=$2
tmpd=/tmp/tmp$$
tmp=$tmpd/tmp.bed
tmpo=$tmpd/tmp.out.bed
tmp2=$tmpd/tmp2.bed
tmpm=$tmpd/tmp.mrg.bed
tmpms=$tmpd/tmp.maxscore.txt

mkdir -p $tmpd

rm -f $tmp
echo "union-ing single-cell sets..."

linenum=1
cat "$fileOfFilenames" | while read peakPath restOfLine
do
    pks=$peakPath
    if [ ! -s "$pks" ]; then
        echo "Error:  File $pks,"
        echo "on line $linenum of $fileOfFilenames, was not found."
        exit 1
    fi

    proj=$(basename "$pks" | cut -f1 -d '.')
    # Each filename must end in .bed.gz or .starch. or .bed.
    is_gz=$(basename "$pks" | grep -c '\.bed\.gz$')
    is_starch=$(basename "$pks" | grep -c '\.starch$')
    is_bed=$(basename "$pks" | grep -c '\.bed$')

    if [ "$is_gz" == "1" ]; then
        zcat "$pks" \
            | awk -v p="$proj" '{print $1"\t"$2"\t"$3"\t"p"\t"$8}' - \
            | bedops -n -1 - "$BLACKLIST" \
            >> "$tmp"
        N=$(zcat "$pks" | wc -l)
    elif [ "$is_starch" == "1" ]; then
        bedops -n -1 "$pks" "$BLACKLIST" \
            | awk -v p="$proj" '{print $1"\t"$2"\t"$3"\t"p"\t"$8}' - \
            >> "$tmp"
        N=$(unstarch "$pks" | wc -l)
    elif [ "$is_bed" == "1" ]; then
        awk -v p="$proj" '{print $1"\t"$2"\t"$3"\t"p"\t"$8}' "$pks" \
            | bedops -n -1 - "$BLACKLIST" \
            >> "$tmp"
        N=$(wc -l < "$pks")
    else
        echo "Error:  Files in $fileOfFilenames must end in .bed.gz, .starch, or .bed."
        echo "File $pks, on line $linenum, does not."
        exit 1
    fi
    printf "%-15s %7d\n" "${proj/-DS.*/}" "$N"
    ((linenum++))
done

echo "sorting..."
sort-bed "$tmp" > "$tmpo"

## Get non-overlapping representatives.
stop=0
while [ "$stop" == 0 ]
do
    echo "merge steps..."
    ## This klugey bit before and after the merging is because we don't want to merge regions that are simply adjacent but not overlapping
    awk '{print $1"\t"$2"\t"$3-1}' $tmpo \
        | bedops -m - \
        | awk '{print $1"\t"$2"\t"$3+1}' - \
        > $tmpm

    bedmap --max $tmpm $tmpo \
        > $tmpms
    awk '{print $1"\t"$2"\t"$3"\tI"}' $tmpm \
        | paste - $tmpms \
        | bedmap --max $tmpo - \
        | paste $tmpo - \
        | awk '($5==$6)' - \
        | cut -f1-4,6 - \
        > $tmp
    ## There may be some ties -- now go through and take one representative from the adjacent ties
    echo "remove adjacent ties..."
    awk '{print $1"\t"$2"\t"$3-1}' $tmp \
        | bedops -m - \
        | awk '{print $1"\t"$2"\t"$3+1}' - \
        | bedmap --multidelim "\t" --echo-map - $tmp \
        | cut -f1-5 \
        > $tmp2
    ## Now add in stuff that doesn't intersect our current set
    bedops -n -0% $tmpo $tmp2 \
        > $tmp
    if [ -s $tmp ]; then
        echo "Adding $(wc -l $tmp | cut -d' ' -f1 -) elements"
        bedops -u $tmp $tmp2 \
            > $tmpo
    else
        cp $tmp2 $tmpo
        stop=1
    fi
done

mv "$tmpo" "$out"
cut -f1-3 "$out" > masterDHSs.bed3
exit
