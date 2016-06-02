#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash run_master_list_simple.sh

ENCODE_AWG_BLACKLIST=wgEncodeHg19ConsensusSignalArtifactRegions.bed

## Make master set of unique DHS positions in the genome.  This will be the union of
## the MCV peaks and all single-cell strict peaks not in MCV zones.
## Also extract the "cell-type predominant" set.  This is everything that's left over from the
## single-cell line union after you subtract off MCV and CTS.
gnom=hg19
out=multi-tissue.master.${gnom}.bed
tmpd=/tmp/tmp$$
tmp=$tmpd/tmp.bed
tmpo=$tmpd/tmp.out.bed
tmp2=$tmpd/tmp2.bed
tmpm=$tmpd/tmp.mrg.bed
tmpms=$tmpd/tmp.maxscore.txt

allpks=(peaks/A549-DS14289.peaks.fdr0.01.hg19.bed.gz
    peaks/AG10803-DS12374.peaks.fdr0.01.hg19.bed.gz
    peaks/AoAF-DS13513.peaks.fdr0.01.hg19.bed.gz
    peaks/CD14-DS17215.peaks.fdr0.01.hg19.bed.gz
    peaks/CD19-DS17186.peaks.fdr0.01.hg19.bed.gz
    peaks/CD20-DS18208.peaks.fdr0.01.hg19.bed.gz
    peaks/CD34-DS12274.peaks.fdr0.01.hg19.bed.gz
    peaks/CD3_CordBlood-DS17706.peaks.fdr0.01.hg19.bed.gz
    peaks/CD3-DS17198.peaks.fdr0.01.hg19.bed.gz
    peaks/CD4-DS17212.peaks.fdr0.01.hg19.bed.gz
    peaks/CD4pos_N-DS14108.peaks.fdr0.01.hg19.bed.gz
    peaks/CD56-DS17189.peaks.fdr0.01.hg19.bed.gz
    peaks/CD8-DS17203.peaks.fdr0.01.hg19.bed.gz
    peaks/fBrain-DS11872.peaks.fdr0.01.hg19.bed.gz
    peaks/fHeart-DS12531.peaks.fdr0.01.hg19.bed.gz
    peaks/fIntestine_Lg-DS17313.peaks.fdr0.01.hg19.bed.gz
    peaks/fKidney-DS20786.peaks.fdr0.01.hg19.bed.gz
    peaks/fLung-DS14724.peaks.fdr0.01.hg19.bed.gz
    peaks/fMuscle_leg-DS20239.peaks.fdr0.01.hg19.bed.gz
    peaks/fPlacenta-DS20346.peaks.fdr0.01.hg19.bed.gz
    peaks/fSkin_fibro_leg_R_quad-DS19943.peaks.fdr0.01.hg19.bed.gz
    peaks/fSpinal_cord-DS20351.peaks.fdr0.01.hg19.bed.gz
    peaks/fStomach-DS17878.peaks.fdr0.01.hg19.bed.gz
    peaks/fThymus-DS20341.peaks.fdr0.01.hg19.bed.gz
    peaks/GM06990-DS7748.peaks.fdr0.01.hg19.bed.gz
    peaks/GM12865-DS12436.peaks.fdr0.01.hg19.bed.gz
    peaks/HAEpiC-DS12663.peaks.fdr0.01.hg19.bed.gz
    peaks/HAh-DS15192.peaks.fdr0.01.hg19.bed.gz
    peaks/HAsp-DS14790.peaks.fdr0.01.hg19.bed.gz
    peaks/HCFaa-DS13480.peaks.fdr0.01.hg19.bed.gz
    peaks/HCF-DS12501.peaks.fdr0.01.hg19.bed.gz
    peaks/HCM-DS12599.peaks.fdr0.01.hg19.bed.gz
    peaks/HCPEpiC-DS12447.peaks.fdr0.01.hg19.bed.gz
    peaks/HEEpiC-DS12763.peaks.fdr0.01.hg19.bed.gz
    peaks/HepG2-DS7764.peaks.fdr0.01.hg19.bed.gz
    peaks/hESCT0-DS11909.peaks.fdr0.01.hg19.bed.gz
    peaks/HFF-DS15115.peaks.fdr0.01.hg19.bed.gz
    peaks/HGF-DS11752.peaks.fdr0.01.hg19.bed.gz
    peaks/HIPEpiC-DS12684.peaks.fdr0.01.hg19.bed.gz
    peaks/HMF-DS13368.peaks.fdr0.01.hg19.bed.gz
    peaks/HMVEC_dBlAd-DS13337.peaks.fdr0.01.hg19.bed.gz
    peaks/HMVEC_dBlNeo-DS13242.peaks.fdr0.01.hg19.bed.gz
    peaks/HMVEC_dLyNeo-DS13150.peaks.fdr0.01.hg19.bed.gz
    peaks/HMVEC_LBl-DS13372.peaks.fdr0.01.hg19.bed.gz
    peaks/HMVEC_LLy-DS13185.peaks.fdr0.01.hg19.bed.gz
    peaks/HPAF-DS13411.peaks.fdr0.01.hg19.bed.gz
    peaks/HPdLF-DS13573.peaks.fdr0.01.hg19.bed.gz
    peaks/HPF-DS13390.peaks.fdr0.01.hg19.bed.gz
    peaks/HRCE-DS10666.peaks.fdr0.01.hg19.bed.gz
    peaks/HSMM-DS14426.peaks.fdr0.01.hg19.bed.gz
    peaks/hTH17-DS11039.peaks.fdr0.01.hg19.bed.gz
    peaks/hTH1-DS7840.peaks.fdr0.01.hg19.bed.gz
    peaks/hTH2-DS17597.peaks.fdr0.01.hg19.bed.gz
    peaks/hTR-DS14702.peaks.fdr0.01.hg19.bed.gz
    peaks/HUVEC-DS10060.peaks.fdr0.01.hg19.bed.gz
    peaks/HVMF-DS13981.peaks.fdr0.01.hg19.bed.gz
    peaks/IMR90-DS13219.peaks.fdr0.01.hg19.bed.gz
    peaks/iPS_19_11-DS15153.peaks.fdr0.01.hg19.bed.gz
    peaks/iTH1-DS18018.peaks.fdr0.01.hg19.bed.gz
    peaks/iTH2-DS17603.peaks.fdr0.01.hg19.bed.gz
    peaks/K562-DS9767.peaks.fdr0.01.hg19.bed.gz
    peaks/LHCN_M2-DS20548.peaks.fdr0.01.hg19.bed.gz
    peaks/M059J-DS20493.peaks.fdr0.01.hg19.bed.gz
    peaks/Mesendoderm-DS19310.peaks.fdr0.01.hg19.bed.gz
    peaks/MSC-DS21042.peaks.fdr0.01.hg19.bed.gz
    peaks/NB4-DS12543.peaks.fdr0.01.hg19.bed.gz
    peaks/NHA-DS12800.peaks.fdr0.01.hg19.bed.gz
    peaks/NHDF_Ad-DS12863.peaks.fdr0.01.hg19.bed.gz
    peaks/NHDF_Neo-DS11923.peaks.fdr0.01.hg19.bed.gz
    peaks/NHLF-DS12829.peaks.fdr0.01.hg19.bed.gz
    peaks/Psoas_Muscle-DS20325.peaks.fdr0.01.hg19.bed.gz
    peaks/RPMI_7951-DS20909.peaks.fdr0.01.hg19.bed.gz
    peaks/SAEC-DS10518.peaks.fdr0.01.hg19.bed.gz
    peaks/Skin_Fibroblasts-DS18224.peaks.fdr0.01.hg19.bed.gz
    peaks/Skin_Keratinocytes-DS18692.peaks.fdr0.01.hg19.bed.gz
    peaks/Skin_Melanocytes-DS18590.peaks.fdr0.01.hg19.bed.gz
    peaks/SkMC-DS11949.peaks.fdr0.01.hg19.bed.gz
    peaks/SKNSH-DS8482.peaks.fdr0.01.hg19.bed.gz
    peaks/Small_Bowel_Mucosa-DS20770.peaks.fdr0.01.hg19.bed.gz
    peaks/T_47D-DS19794.peaks.fdr0.01.hg19.bed.gz
    peaks/Trophoblast-DS19317.peaks.fdr0.01.hg19.bed.gz
    peaks/vHMEC-DS18406.peaks.fdr0.01.hg19.bed.gz)

mkdir -p "$tmpd"

rm -f "$tmp"
echo "union-ing single-cell sets..."

for pks in ${allpks[*]}
do
    proj=$(basename "$pks" | cut -f1 -d '.')
    zcat "$pks" \
    | awk -v p="$proj" '{print $1"\t"$2"\t"$3"\t"p"\t"$8}' - \
    | bedops -n -1 - "$ENCODE_AWG_BLACKLIST" \
    >> "$tmp"
    N=$(zcat "$pks" | wc -l)
    printf "%-15s %7d\n" "${proj/-DS.*/}" "$N"
done

echo "sorting..."
sort-bed $tmp > $tmpo

## Get non-overlapping representatives.
stop=0
while [ $stop == 0 ]
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
    if [ -s $tmp ]
    then
        echo "Adding $(wc -l $tmp | cut -d' ' -f1 -) elements"
    bedops -u $tmp $tmp2 \
        > $tmpo
    else
    cp $tmp2 $tmpo
    stop=1
    fi
done

mv $tmpo $out
exit
