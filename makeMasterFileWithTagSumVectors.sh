#! /bin/bash
# qsub -cwd -N JOB_ID -S /bin/bash makeMasterFileWithTagSumVectors.sh

source ~/.bashrc

PID=$$
TEMP1=tempfile1_${PID}.bed
TEMP2=tempfile2_${PID}.bed
outfile=masterDHSsAndTagCounts_82_hg19.bed4

TT=totalDeepSeqTagCountsPerCelltype.txt

files=(tagDensities/A549-DS14289_tagDensityInMasterDHSs.bed4
tagDensities/AG10803-DS12374_tagDensityInMasterDHSs.bed4
tagDensities/AoAF-DS13513_tagDensityInMasterDHSs.bed4
tagDensities/CD14-DS17215_tagDensityInMasterDHSs.bed4
tagDensities/CD19-DS17186_tagDensityInMasterDHSs.bed4
tagDensities/CD20-DS18208_tagDensityInMasterDHSs.bed4
tagDensities/CD34-DS12274_tagDensityInMasterDHSs.bed4
tagDensities/CD3_CordBlood-DS17706_tagDensityInMasterDHSs.bed4
tagDensities/CD3-DS17198_tagDensityInMasterDHSs.bed4
tagDensities/CD4-DS17212_tagDensityInMasterDHSs.bed4
tagDensities/CD4pos_N-DS14108_tagDensityInMasterDHSs.bed4
tagDensities/CD56-DS17189_tagDensityInMasterDHSs.bed4
tagDensities/CD8-DS17203_tagDensityInMasterDHSs.bed4
tagDensities/fBrain-DS11872_tagDensityInMasterDHSs.bed4
tagDensities/fHeart-DS12531_tagDensityInMasterDHSs.bed4
tagDensities/fIntestine_Lg-DS17313_tagDensityInMasterDHSs.bed4
tagDensities/fKidney-DS20786_tagDensityInMasterDHSs.bed4
tagDensities/fLung-DS14724_tagDensityInMasterDHSs.bed4
tagDensities/fMuscle_leg-DS20239_tagDensityInMasterDHSs.bed4
tagDensities/fPlacenta-DS20346_tagDensityInMasterDHSs.bed4
tagDensities/fSkin_fibro_leg_R_quad-DS19943_tagDensityInMasterDHSs.bed4
tagDensities/fSpinal_cord-DS20351_tagDensityInMasterDHSs.bed4
tagDensities/fStomach-DS17878_tagDensityInMasterDHSs.bed4
tagDensities/fThymus-DS20341_tagDensityInMasterDHSs.bed4
tagDensities/GM06990-DS7748_tagDensityInMasterDHSs.bed4
tagDensities/GM12865-DS12436_tagDensityInMasterDHSs.bed4
tagDensities/HAEpiC-DS12663_tagDensityInMasterDHSs.bed4
tagDensities/HAh-DS15192_tagDensityInMasterDHSs.bed4
tagDensities/HAsp-DS14790_tagDensityInMasterDHSs.bed4
tagDensities/HCFaa-DS13480_tagDensityInMasterDHSs.bed4
tagDensities/HCF-DS12501_tagDensityInMasterDHSs.bed4
tagDensities/HCM-DS12599_tagDensityInMasterDHSs.bed4
tagDensities/HCPEpiC-DS12447_tagDensityInMasterDHSs.bed4
tagDensities/HEEpiC-DS12763_tagDensityInMasterDHSs.bed4
tagDensities/HepG2-DS7764_tagDensityInMasterDHSs.bed4
tagDensities/hESCT0-DS11909_tagDensityInMasterDHSs.bed4
tagDensities/HFF-DS15115_tagDensityInMasterDHSs.bed4
tagDensities/HGF-DS11752_tagDensityInMasterDHSs.bed4
tagDensities/HIPEpiC-DS12684_tagDensityInMasterDHSs.bed4
tagDensities/HMF-DS13368_tagDensityInMasterDHSs.bed4
tagDensities/HMVEC_dBlAd-DS13337_tagDensityInMasterDHSs.bed4
tagDensities/HMVEC_dBlNeo-DS13242_tagDensityInMasterDHSs.bed4
tagDensities/HMVEC_dLyNeo-DS13150_tagDensityInMasterDHSs.bed4
tagDensities/HMVEC_LBl-DS13372_tagDensityInMasterDHSs.bed4
tagDensities/HMVEC_LLy-DS13185_tagDensityInMasterDHSs.bed4
tagDensities/HPAF-DS13411_tagDensityInMasterDHSs.bed4
tagDensities/HPdLF-DS13573_tagDensityInMasterDHSs.bed4
tagDensities/HPF-DS13390_tagDensityInMasterDHSs.bed4
tagDensities/HRCE-DS10666_tagDensityInMasterDHSs.bed4
tagDensities/HSMM-DS14426_tagDensityInMasterDHSs.bed4
tagDensities/hTH17-DS11039_tagDensityInMasterDHSs.bed4
tagDensities/hTH1-DS7840_tagDensityInMasterDHSs.bed4
tagDensities/hTH2-DS17597_tagDensityInMasterDHSs.bed4
tagDensities/hTR-DS14702_tagDensityInMasterDHSs.bed4
tagDensities/HUVEC-DS10060_tagDensityInMasterDHSs.bed4
tagDensities/HVMF-DS13981_tagDensityInMasterDHSs.bed4
tagDensities/IMR90-DS13219_tagDensityInMasterDHSs.bed4
tagDensities/iPS_19_11-DS15153_tagDensityInMasterDHSs.bed4
tagDensities/iTH1-DS18018_tagDensityInMasterDHSs.bed4
tagDensities/iTH2-DS17603_tagDensityInMasterDHSs.bed4
tagDensities/K562-DS9767_tagDensityInMasterDHSs.bed4
tagDensities/LHCN_M2-DS20548_tagDensityInMasterDHSs.bed4
tagDensities/M059J-DS20493_tagDensityInMasterDHSs.bed4
tagDensities/Mesendoderm-DS19310_tagDensityInMasterDHSs.bed4
tagDensities/MSC-DS21042_tagDensityInMasterDHSs.bed4
tagDensities/NB4-DS12543_tagDensityInMasterDHSs.bed4
tagDensities/NHA-DS12800_tagDensityInMasterDHSs.bed4
tagDensities/NHDF_Ad-DS12863_tagDensityInMasterDHSs.bed4
tagDensities/NHDF_Neo-DS11923_tagDensityInMasterDHSs.bed4
tagDensities/NHLF-DS12829_tagDensityInMasterDHSs.bed4
tagDensities/Psoas_Muscle-DS20325_tagDensityInMasterDHSs.bed4
tagDensities/RPMI_7951-DS20909_tagDensityInMasterDHSs.bed4
tagDensities/SAEC-DS10518_tagDensityInMasterDHSs.bed4
tagDensities/Skin_Fibroblasts-DS18224_tagDensityInMasterDHSs.bed4
tagDensities/Skin_Keratinocytes-DS18692_tagDensityInMasterDHSs.bed4
tagDensities/Skin_Melanocytes-DS18590_tagDensityInMasterDHSs.bed4
tagDensities/SkMC-DS11949_tagDensityInMasterDHSs.bed4
tagDensities/SKNSH-DS8482_tagDensityInMasterDHSs.bed4
tagDensities/Small_Bowel_Mucosa-DS20770_tagDensityInMasterDHSs.bed4
tagDensities/T_47D-DS19794_tagDensityInMasterDHSs.bed4
tagDensities/Trophoblast-DS19317_tagDensityInMasterDHSs.bed4
tagDensities/vHMEC-DS18406_tagDensityInMasterDHSs.bed4)

MAX_COUNT=0
for file in ${files[*]}
do
    celltype=`basename $file _tagDensityInMasterDHSs.bed4`
    count=`grep -w $celltype $TT | cut -f2`
    if [ "$count" -gt "$MAX_COUNT" ]; then
	MAX_COUNT=$count
    fi
done
if [ "$MAX_COUNT" == "0" ]; then
    echo -e "Error:  Failed to find maximum genomewide tag count. Exiting...."
    exit
fi

file=A549-DS14289_tagDensityInMasterDHSs.bed4
celltype=`basename $file _tagDensityInMasterDHSs.bed4`
count=`grep -w $celltype $TT | cut -f2`
if [ "$count" == "" ]; then
    echo -e "Failed to find $celltype in $TT."
    exit
fi

awk -v count=$count -v maxCount=$MAX_COUNT 'BEGIN{OFS="\t"}{print $1,$2,$3,int(($4*1.0*maxCount/count) + 0.5);}' $file > $TEMP2

for file in ${files[*]}
do
    if [ "$file" == "A549-DS14289_tagDensityInMasterDHSs.bed4" ]; then
	continue
    fi
    celltype=`basename $file _tagDensityInMasterDHSs.bed4`
    count=`grep -w $celltype $TT | cut -f2`
    if [ "$count" == "" ]; then
	echo -e "Failed to find $celltype in $TT."
	exit
    fi
    awk -v count=$count -v maxCount=$MAX_COUNT 'BEGIN{OFS="\t"}{print int(($4*1.0*maxCount/count) + 0.5);}' $file \
	| paste -d "," $TEMP2 - \
	> $TEMP1
    mv $TEMP1 $TEMP2
done
mv $TEMP2 $outfile
rm -f $TEMP1
exit
