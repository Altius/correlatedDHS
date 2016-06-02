#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash getTagDensitiesInMasterListDHSs_pseudoParallel.sh min_idx max_idx

if [ "$1" == "" ] || [ "$2" == "" ]; then
    echo -e "Usage:  $0 minIdxInclusive maxIdxInclusive"
    exit
fi

ENCODE_AWG_BLACKLIST=wgEncodeHg19ConsensusSignalArtifactRegions.bed

MIN_IDX=$1
MAX_IDX=$2

OUTDIR=tagDensities
mkdir -p $OUTDIR

EXE=gchr
MASTER_DHSs=masterDHSs.bed3
files=(density/A549-DS14289.75.20.uniques-density.36.hg19.bed.jarch
    density/AG10803-DS12374.75.20.uniques-density.36.hg19.bed.jarch
    density/AoAF-DS13513.75.20.uniques-density.36.hg19.bed.jarch
    density/CD14-DS17215.75.20.uniques-density.36.hg19.bed.jarch
    density/CD19-DS17186.75.20.uniques-density.36.hg19.bed.jarch
    density/CD20-DS18208.75.20.uniques-density.36.hg19.bed.jarch
    density/CD34-DS12274.75.20.uniques-density.36.hg19.bed.jarch
    density/CD3_CordBlood-DS17706.75.20.uniques-density.36.hg19.bed.jarch
    density/CD3-DS17198.75.20.uniques-density.36.hg19.bed.jarch
    density/CD4-DS17212.75.20.uniques-density.36.hg19.bed.jarch
    density/CD4pos_N-DS14108.75.20.uniques-density.36.hg19.bed.jarch
    density/CD56-DS17189.75.20.uniques-density.36.hg19.bed.jarch
    density/CD8-DS17203.75.20.uniques-density.36.hg19.bed.jarch
    density/fBrain-DS11872.75.20.uniques-density.36.hg19.bed.jarch
    density/fHeart-DS12531.75.20.uniques-density.36.hg19.bed.jarch
    density/fIntestine_Lg-DS17313.75.20.uniques-density.36.hg19.bed.jarch
    density/fKidney-DS20786.75.20.uniques-density.36.hg19.bed.jarch
    density/fLung-DS14724.75.20.uniques-density.36.hg19.bed.jarch
    density/fMuscle_leg-DS20239.75.20.uniques-density.36.hg19.bed.jarch
    density/fPlacenta-DS20346.75.20.uniques-density.36.hg19.bed.jarch
    density/fSkin_fibro_leg_R_quad-DS19943.75.20.uniques-density.36.hg19.bed.jarch
    density/fSpinal_cord-DS20351.75.20.uniques-density.36.hg19.bed.jarch
    density/fStomach-DS17878.75.20.uniques-density.36.hg19.bed.jarch
    density/fThymus-DS20341.75.20.uniques-density.36.hg19.bed.jarch
    density/GM06990-DS7748.75.20.uniques-density.27.hg19.bed.jarch
    density/GM12865-DS12436.75.20.uniques-density.36.hg19.bed.jarch
    density/HAEpiC-DS12663.75.20.uniques-density.36.hg19.bed.jarch
    density/HAh-DS15192.75.20.uniques-density.36.hg19.bed.jarch
    density/HAsp-DS14790.75.20.uniques-density.36.hg19.bed.jarch
    density/HCFaa-DS13480.75.20.uniques-density.36.hg19.bed.jarch
    density/HCF-DS12501.75.20.uniques-density.36.hg19.bed.jarch
    density/HCM-DS12599.75.20.uniques-density.36.hg19.bed.jarch
    density/HCPEpiC-DS12447.75.20.uniques-density.36.hg19.bed.jarch
    density/HEEpiC-DS12763.75.20.uniques-density.36.hg19.bed.jarch
    density/HepG2-DS7764.75.20.uniques-density.27.hg19.bed.jarch
    density/hESCT0-DS11909.75.20.uniques-density.36.hg19.bed.jarch
    density/HFF-DS15115.75.20.uniques-density.36.hg19.bed.jarch
    density/HGF-DS11752.75.20.uniques-density.36.hg19.bed.jarch
    density/HIPEpiC-DS12684.75.20.uniques-density.36.hg19.bed.jarch
    density/HMF-DS13368.75.20.uniques-density.36.hg19.bed.jarch
    density/HMVEC_dBlAd-DS13337.75.20.uniques-density.36.hg19.bed.jarch
    density/HMVEC_dBlNeo-DS13242.75.20.uniques-density.36.hg19.bed.jarch
    density/HMVEC_dLyNeo-DS13150.75.20.uniques-density.36.hg19.bed.jarch
    density/HMVEC_LBl-DS13372.75.20.uniques-density.36.hg19.bed.jarch
    density/HMVEC_LLy-DS13185.75.20.uniques-density.36.hg19.bed.jarch
    density/HPAF-DS13411.75.20.uniques-density.36.hg19.bed.jarch
    density/HPdLF-DS13573.75.20.uniques-density.36.hg19.bed.jarch
    density/HPF-DS13390.75.20.uniques-density.36.hg19.bed.jarch
    density/HRCE-DS10666.75.20.uniques-density.36.hg19.bed.jarch
    density/HSMM-DS14426.75.20.uniques-density.36.hg19.bed.jarch
    density/hTH17-DS11039.75.20.uniques-density.36.hg19.bed.jarch
    density/hTH1-DS7840.75.20.uniques-density.27.hg19.bed.jarch
    density/hTH2-DS17597.75.20.uniques-density.36.hg19.bed.jarch
    density/hTR-DS14702.75.20.uniques-density.36.hg19.bed.jarch
    density/HUVEC-DS10060.75.20.uniques-density.36.hg19.bed.jarch
    density/HVMF-DS13981.75.20.uniques-density.36.hg19.bed.jarch
    density/IMR90-DS13219.75.20.uniques-density.36.hg19.bed.jarch
    density/iPS_19_11-DS15153.75.20.uniques-density.36.hg19.bed.jarch
    density/iTH1-DS18018.75.20.uniques-density.36.hg19.bed.jarch
    density/iTH2-DS17603.75.20.uniques-density.36.hg19.bed.jarch
    density/K562-DS9767.75.20.uniques-density.27.hg19.bed.jarch
    density/LHCN_M2-DS20548.75.20.uniques-density.36.hg19.bed.jarch
    density/M059J-DS20493.75.20.uniques-density.36.hg19.bed.jarch
    density/Mesendoderm-DS19310.75.20.uniques-density.36.hg19.bed.jarch
    density/MSC-DS21042.75.20.uniques-density.36.hg19.bed.jarch
    density/NB4-DS12543.75.20.uniques-density.36.hg19.bed.jarch
    density/NHA-DS12800.75.20.uniques-density.36.hg19.bed.jarch
    density/NHDF_Ad-DS12863.75.20.uniques-density.36.hg19.bed.jarch
    density/NHDF_Neo-DS11923.75.20.uniques-density.36.hg19.bed.jarch
    density/NHLF-DS12829.75.20.uniques-density.36.hg19.bed.jarch
    density/Psoas_Muscle-DS20325.75.20.uniques-density.36.hg19.bed.jarch
    density/RPMI_7951-DS20909.75.20.uniques-density.36.hg19.bed.jarch
    density/SAEC-DS10518.75.20.uniques-density.36.hg19.bed.jarch
    density/Skin_Fibroblasts-DS18224.75.20.uniques-density.36.hg19.bed.jarch
    density/Skin_Keratinocytes-DS18692.75.20.uniques-density.36.hg19.bed.jarch
    density/Skin_Melanocytes-DS18590.75.20.uniques-density.36.hg19.bed.jarch
    density/SkMC-DS11949.75.20.uniques-density.36.hg19.bed.jarch
    density/SKNSH-DS8482.75.20.uniques-density.27.hg19.bed.jarch
    density/Small_Bowel_Mucosa-DS20770.75.20.uniques-density.36.hg19.bed.jarch
    density/T_47D-DS19794.75.20.uniques-density.36.hg19.bed.jarch
    density/Trophoblast-DS19317.75.20.uniques-density.36.hg19.bed.jarch
    density/vHMEC-DS18406.75.20.uniques-density.36.hg19.bed.jarch)

i=$MIN_IDX
while [ "$i" -le "$MAX_IDX" ]; do
# for file in ${files[*]}
# do
    file=${files[i]}
    fnameStub=$(basename "$file" | cut -f1 -d '.')
    outfile=$OUTDIR/${fnameStub}_tagDensityInMasterDHSs.bed4
    # No call to sort-bed needed after the awk statement in this case.
    # All fractions will be a multiple of 0.25 in this case, hence --prec 2.
    #    $EXE $file \
    #        | bedmap --echo --bases-uniq-f --echo-map - $MASTER_DHSs \
    #        | awk -F "|" '{if($2!=0){split($1,x,"\t");print $3"\tid\t"$2*x[5];}}' \
    #        | bedmap --delim "\t" --prec 2 --echo --sum $MASTER_DHSs - \
    #        > $outfile

    # Whoops!  There are many DHSs that are adjacent to each other,
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

    "$EXE" "$file" \
        | bedops -n -1 - "$ENCODE_AWG_BLACKLIST" \
        | bedmap --count --echo --bases-uniq-f --echo-map - "$MASTER_DHSs" \
        | awk -F "|" \
            '{
                 split($2, x, "\t");
                 dens = x[5] / 7.5;
                 if(1==$1) {
                     print $4 "\tid\t" dens * $3;
                 } else if (2==$1) {
                     split($4, y, ";");
                     split(y[1], z, "\t");
                     print y[1] "\tid\t" dens * (z[3] - x[2]) * 0.05;
                     split(y[2], z, "\t");
                     print y[2] "\tid\t" dens * (x[3] - z[2]) * 0.05;
                 }
             }' \
        | bedmap --delim "\t" --prec 3 --echo --sum "$MASTER_DHSs" - \
        > "$outfile"
    ((i++))
done
exit

