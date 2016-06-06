#! /bin/bash
# doEverything.sh

NUMFILES=${NUMFILES:-82}

submit(){
    qsub -cwd -V -S /bin/bash "$@"
}

set -e -o pipefail # quit when any section of any pipe, or single command, fails

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  $0 fileOfPeakfileNames.txt fileOfDensityfileNames.txt"
    exit 1
fi

peakfileNames=$1
densityfileNames=$2

#qsub -cwd -N RunML -S run_master_list_simple__receiveInputFileOfFilenames.sh $peakfileNames
# Hmmm, need to check exit status....

PATH=$(dirname $0)/$PATH

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

    run_master_list_simple.sh
    cut -f1-3 multi-tissue.master.hg19.bed > masterDHSs.bed3

    date
__SETUP__

submit -hold_jid "$SETUP_JOB" -N "$DENS_JOB" -t "$NUMFILES" <<__DENSITY__
    set -x
    date

    getTagDensitiesInMasterListDHSs_pseudoParallel.sh "$SGE_TASK_ID" "$SGE_TASK_ID"

    date
__DENSITY__

submit -hold_jid "$DENS_JOB" -N "$FINISH_JOB" <<__FINISH__
    set -x
    date
    makeMasterFileWithTagSumVectors.sh

    date
    getPromoterDHSs \
        EH_v2_TxStarts.bed6 masterDHSsAndTagCounts_82_hg19.bed4 \
        10000 2500 20 100 \
        promOutfile.bed13

    date
    awk 'BEGIN{OFS="\t"}{if($13!="NA"){print $5,$6,$7,$4,".",$8,$1,$2,$3;}}'
    promOutfile.bed13 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9

    date
    cut -f1-4 promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsWithGeneNames_EHv2_82celltypes_hg19.bed4

    date
    bedmap --exact --delim "\t" --echo --echo-map-id \
        promDHSsWithGeneNames_EHv2_82celltypes_hg19.bed4 \
        masterDHSsAndTagCounts_82_hg19.bed4 \
        > promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5

    date
    cut -f1-3 promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5 \
        | uniq \
        | bedops -n -1 masterDHSsAndTagCounts_82_hg19.bed4 - \
        > nonprom_temp.bed4

    date
    closest-features --shortest --dist nonprom_temp.bed4 \
        promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9 \
        | awk -F "|" '{dist=$3;if(dist<0){dist=-dist;}if(dist>=1000){print $1;}}' \
        > nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS_82celltypes.bed4

    date
    bedmap --range 500000 --skip-unmapped --echo --echo-map \
        promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5 \
        nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS_82celltypes.bed4 \
        > inputForCorrelationCalculations_500kb_82celltypes.bed

    date
    calcCorrelations.sh \
        inputForCorrelationCalculations_500kb_82celltypes.bed \
        0.7 \
        corrs_promsFirst_EHv2_all_82celltypes_500kb.bed8 \
        corrs_promsFirst_EHv2_above0.7_82celltypes_500kb.bed8 \
        corrs_distalsFirst_EHv2_all_82celltypes_500kb.bed8 \
        corrs_distalsFirst_EHv2_above0.7_82celltypes_500kb.bed8

    date
    rm inputForCorrelationCalculations_500kb_82celltypes.bed
    __FINISH__
