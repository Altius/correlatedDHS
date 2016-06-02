#! /bin/bash
# qsub -cwd -N GetProms -S /bin/bash getPromoterDHSs_TryToIncreaseQuality.sh

source ~/.bashrc

if [ "$1" == "" ] || [ "$2" == "" ] || [ "$3" == "" ] || [ "$4" == "" ] || [ "$5" == "" ] || [ "$6" == "" ] || [ "$7" == "" ]; then
    echo -e "Usage: "$0" TSSfile.bed6 DHSsAndTheirTagCounts.bed4 UpstreamThresholdInBp DownstreamThresholdInBp minTagCount tiebreakerTagCount Outfile.bed13"
    echo -e "Generally, the nearest DHS in the 5' direction of the TSS will get called as \"the promoter DHS.\""
    echo -e "If this DHS is too far upstream from the TSS (>UpstreamThresholdInBp, e.g. 10000) or too weak (minTagCount (e.g. 20) not attained in any cell type),"
    echo -e "then that DHS won't get called as the promoter DHS."
    echo -e "When this happens, the nearest DHS in the 3' direction is tested."
    echo -e "If this DHS is too far from the TSS (>DownstreamThresholdInBp, e.g. 1500) or too weak (minTagCount not attained in any cell type),"
    echo -e "then no promoter DHS will get called for the TSS (\"NA\")."
    echo -e "If the nearest DHS in the 3' direction is nearer to the TSS than the nearest DHS in the 5' direction,"
    echo -e "then if it's much nearer (at least 3kb nearer) and is strong enough, it \"wins the tiebreaker\" and gets called as the promoter DHS."
    echo -e "If both DHSs are relatively near the TSS and the 3' one is very strong (>tiebreakerTagCount, e.g. 100) while the 5' one is not (minTagCount < count < tiebreakerCount),"
    echo -e "then the 3' DHS wins the tiebreaker.  Otherwise the 5' DHS gets called."
    echo -e "The final field on each line specifies which DHS was called, when applicable."
    exit
fi

PID=$$
GENE_TSS_FILE=$1
DHS_INPUT_FILE=$2
UPSTREAM_THRESHOLD_IN_BP=$3
DOWNSTREAM_THRESHOLD_IN_BP=$4
MIN_TAGCOUNT=$5
TIEBREAKER_TAGCOUNT=$6
OUTPUT_FILE=$7

if [ ! -f $GENE_TSS_FILE ]; then
    echo -e "Unable to find file "$GENE_TSS_FILE"."
    exit
fi

if [ ! -f $DHS_INPUT_FILE ]; then
    echo -e "Unable to find file "$DHS_INPUT_FILE"."
    exit
fi

# $GENE_TSS_FILE must have at least 6 fields, with the strand in field 6.
if [ ! `head -n 1 $GENE_TSS_FILE | cut -f6` ]; then
    echo -e "File "$GENE_TSS_FILE" must have the strand in field 6."
    exit
fi

# DHS_INPUT_FILE must have exactly 4 fields, with the comma-delimited tag counts in field 4.
if [ ! `head -n 1 $DHS_INPUT_FILE | cut -f4` ] || [ `head -n 1 $DHS_INPUT_FILE | cut -f5` ]; then
    echo -e "File "$DHS_INPUT_FILE" must have exactly 4 columns, with the comma-delimited tag counts in field 4."
    exit
fi

# Just in case there are >6 fields in the input TSS file,
# start by calling cut -f1-6.

echo -e "Executing: $0 $GENE_TSS_FILE $DHS_INPUT_FILE $UPSTREAM_THRESHOLD_IN_BP $DOWNSTREAM_THRESHOLD_IN_BP $MIN_TAGCOUNT $TIEBREAKER_TAGCOUNT $OUTPUT_FILE"

# Note:  When the distance between a TSS and the closest DHS is 0,
# closest-features always reports that DHS as the "left" DHS.
# So, regardless of the TSS's strand, whenever the left DHS has distance=0,
# we need to select it instead of the right DHS.
# See the usage statement at the top of this script for tiebreaking scenarios.
cut -f1-6 $GENE_TSS_FILE | closest-features --delim '\t' --dist - $DHS_INPUT_FILE \
	| gawk -v upThresh=$UPSTREAM_THRESHOLD_IN_BP -v downThresh=$DOWNSTREAM_THRESHOLD_IN_BP -v minTagCount=$MIN_TAGCOUNT -v tiebreakerTagCount=$TIEBREAKER_TAGCOUNT 'BEGIN{FS="\t"; OFS="\t"}{ \
      if ($8=="NA") { \
          Ldist=999999999; \
          maxLcount=0; \
          if ($10=="NA") { \
             Rdist=999999999; \
             maxRcount=0;} \
          else { \
             Rdist=$13; \
             Rchr=$9; \
             Rbeg=$10; \
             Rend=$11; \
             split($12,x,","); \
             asort(x); \
             maxRcount=x[length(x)];} \
      } \
      else { \
         Ldist = -$11; \
         Lchr=$7; \
         Lbeg=$8; \
         Lend=$9; \
         split($10,x,","); \
         asort(x); \
         maxLcount=x[length(x)]; \
         if ($13=="NA") { \
            Rdist=999999999; \
            maxRcount=0;} \
         else { \
             Rdist=$16; \
             Rchr=$12; \
             Rbeg=$13; \
             Rend=$14; \
             split($15,x,","); \
             asort(x); \
             maxRcount=x[length(x)];} \
      } \
      if ($6=="+" || Ldist==0) { \
         dist=Ldist; \
         chr=Lchr; \
         beg=Lbeg; \
         end=Lend; \
         maxCount=maxLcount; \
         altDist=Rdist; \
         altChr=Rchr; \
         altBeg=Rbeg; \
         altEnd=Rend; \
         altMaxCount=maxRcount;} \
      else { \
         dist=Rdist; \
         chr=Rchr; \
         beg=Rbeg; \
         end=Rend; \
         maxCount=maxRcount; \
         altDist=Ldist; \
         altChr=Lchr; \
         altBeg=Lbeg; \
         altEnd=Lend; \
         altMaxCount=maxLcount; \
      } \
      if (dist<=upThresh && maxCount >= minTagCount) { \
         if (altDist<=downThresh && altMaxCount >= minTagCount && altDist < dist && \
             ((altDist + 3000 < dist) || (altMaxCount >= tiebreakerThresh && maxCount < tiebreakerThresh))) { \
            print $1,$2,$3,$4,altChr,altBeg,altEnd,$6,dist,altDist,maxCount,altMaxCount,"3prime"; } \
         else { \
            print $1,$2,$3,$4,chr,beg,end,$6,dist,altDist,maxCount,altMaxCount,"5prime"; } \
      } \
      else { \
         if (altDist<=downThresh && altMaxCount >= minTagCount) { \
            print $1,$2,$3,$4,altChr,altBeg,altEnd,$6,dist,altDist,maxCount,altMaxCount,"3prime"; } \
         else { \
            print $1,$2,$3,$4,"NA","NA","NA",$6,dist,altDist,maxCount,altMaxCount,"NA";} \
      } \
   }' > $OUTPUT_FILE

exit
