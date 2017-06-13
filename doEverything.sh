#! /bin/bash
# doEverything.sh

submit(){
    qsub -cwd -V -S /bin/bash "$@"
}

set -e -o pipefail # quit when any section of any pipe, or single command, fails

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  $0 fileOfPeakDensityAndCounts.txt TSSfile.bed6" "[blacklistfile]"
    echo -e "\tfileOfPeakDensityAndCounts.txt should contain one sample per line in the following format:"
    echo -e "\tColumn 1: /path/to/peak/file"
    echo -e "\tColumn 2: /path/to/density/file"
    echo -e "\tColumn 3: Total mapped tags"
    echo -e "\tBlacklist file is optional"
    exit 1
fi

sampleFile=$1
tssFile=$2

if [ ! -f "$sampleFile" ]; then
    echo -e "Unable to find file \"$sampleFile\"."
    exit
fi

if [ ! -f "$tssFile" ]; then
    echo -e "Unable to find file \"$tssFile\"."
    exit
fi
#Check for optional blacklist file, set blacklist to /dev/null if not supplied
if [ ! -z "$3" ]; then
    blacklist=$3
    if [ ! -f "$blacklist" ]; then
        echo -e "Unable to find file \"$blacklist\"."
        exit
    fi
else
    blacklist=/dev/null
fi

#Get the number of samples
NUMFILES=$(cat $sampleFile | wc -l)

PID=$$
SAMPLE_NAMES=sampleNames_${PID}.txt
RAW_COUNTS=rawCounts_${PID}.txt

cut -f3 $sampleFile > $RAW_COUNTS
cat "$sampleFile" | while read peakPath restOfLine
do
    basename "$peakPath" | cut -f1 -d '.' >> $SAMPLE_NAMES
done
paste $SAMPLE_NAMES $RAW_COUNTS > totalTagCountsPerSample.txt
rm -f $RAW_COUNTS $SAMPLE_NAMES

#qsub -cwd -N RunML -S run_master_list_simple__receiveInputFileOfFilenames.sh $peakfileNames
# Hmmm, need to check exit status....

PATH=$(dirname $0):$PATH
PATH=$(dirname $0)/bin:$PATH

if ! which bedops 2>&1 >/dev/null ; then
    echo "Bedops must be installed and available on the PATH" >&2
    exit 2
fi

SETUP_JOB=".setupdhs"
DENS_JOB=".densitydhs"
FINISH_JOB=".finishdhs"

submit -N "$SETUP_JOB" <<__SETUP__
    set -x
    date

    run_master_list_simple__receiveInputFileOfFilenames.sh "$sampleFile" multi-tissue.master.bed "$blacklist"

    date
__SETUP__


submit -hold_jid "$SETUP_JOB" -N "$DENS_JOB" <<__DENSITY__
    set -x
    date

    getTagDensitiesInMasterListDHSs_pseudoParallel_receiveFileOfFilenames.sh $sampleFile $blacklist

    date
__DENSITY__


submit -hold_jid "$DENS_JOB" -N "$FINISH_JOB" <<__FINISH__
    set -x
    date
    makeMasterFileWithTagSumVectors.sh

    date
    getPromoterDHSs.sh \
        $tssFile masterDHSsAndTagCounts.bed4 \
        10000 2500 20 100 \
        promOutfile.bed13

    date
    awk 'BEGIN{OFS="\t"}{if(\$13!="NA"){print \$5,\$6,\$7,\$4,".",\$8,\$1,\$2,\$3;}}' \
    promOutfile.bed13 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsAndTheirTSSs.bed9

    date
    cut -f1-4 promDHSsAndTheirTSSs.bed9 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsWithGeneNames.bed4

    date
    bedmap --exact --delim "\t" --echo --echo-map-id \
        promDHSsWithGeneNames.bed4 \
        masterDHSsAndTagCounts.bed4 \
        > promDHSsWithGeneNamesAndTagCounts.bed5

    date
    cut -f1-3 promDHSsWithGeneNamesAndTagCounts.bed5 \
        | uniq \
        | bedops -n -1 masterDHSsAndTagCounts.bed4 - \
        > nonprom_temp.bed4

    date
    closest-features --shortest --dist nonprom_temp.bed4 \
        promDHSsAndTheirTSSs.bed9 \
        | awk -F "|" '{dist=\$3;if(dist<0){dist=-dist;}if(dist>=1000){print \$1;}}' \
        > nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS.bed4

    date
    bedmap --range 500000 --skip-unmapped --echo --echo-map \
        promDHSsWithGeneNamesAndTagCounts.bed5 \
        nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS.bed4 \
        > inputForCorrelationCalculations_500kb.bed

    date
    calcCorrelations.sh \
        inputForCorrelationCalculations_500kb.bed \
        0.7 \
        corrs_promsFirst_all_celltypes_500kb.bed8 \
        corrs_promsFirst_above0.7_celltypes_500kb.bed8 \
        corrs_distalsFirst_all_celltypes_500kb.bed8 \
        corrs_distalsFirst_above0.7_celltypes_500kb.bed8

    date
__FINISH__
