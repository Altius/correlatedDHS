# Description

This is a set of scripts and an executable to calculate DHS correlation between samples or cell-types.

# How to Use

## Preliminary Steps.

If you have not already installed BEDOPS, you'll need to do so, from
<https://github.com/bedops/bedops>

Compile the code with:

    make

## Easy usage:

If you have an SGE cluster available (`qsub`), run `./doEverything.sh` and wait for the jobs to finish.


# How it works, step-by-step

The `./doEverything.sh` script does the following steps for you - if you want to see the reasoning for any steps, look here:

## Creating the "master list" file.

Start by running this script (preferably on a node of a computer cluster):

    run_master_list_simple.sh

After that has completed and created the file multi-tissue.master.hg19.bed, run
the following command from the command line:

    cut -f1-3 multi-tissue.master.hg19.bed > masterDHSs.bed3

Then run the following script, which gets the DNase signal strength in each cell
type at each DHS. This script requires two integer parameters, which are indexes
into the array of cell types (zero-based). This allows you to process multiple
cell types simultaneously on a computer cluster.

    getTagDensitiesInMasterListDHSs_pseudoParallel.sh idxMin idxMax

For example, you could submit the following jobs to a cluster:

    qsub -cwd -N Get00to07 -S /bin/bash getTagDensitiesInMasterListDHSs_pseudoParallel.sh 0 7
    qsub -cwd -N Get08to15 -S /bin/bash getTagDensitiesInMasterListDHSs_pseudoParallel.sh 8 15
    qsub -cwd -N Get16to23 -S /bin/bash getTagDensitiesInMasterListDHSs_pseudoParallel.sh 16 23

etc. etc. The first job would process the first 8 cell types (0-7), the next one
would process the next 8 cell types (8-15), etc., simultaneously, with one job
per node of a computer cluster. There are 82 cell types total (0-81).

IMPORTANT NOTE:  The script `getTagDensitiesInMasterListDHSs_pseudoParallel.sh`
calls the script gchr (included), which is a "C shell" script. If
`getTagDensitiesInMasterListDHSs_pseudoParallel.sh` fails to execute gchr, ask
Michael for guidance on how to set things up. (All scripts in this directory are
designed to run in Bourne shell; they might run fine in other shells too; I
don't know offhand what issues might arise.)

After all of those jobs have completed and produced 82 files in subdirectory
tagDensities, run the following script:

    makeMasterFileWithTagSumVectors.sh

This will create the file `masterDHSsAndTagCounts_82_hg19.bed4`, in which fields
1-3 are the coordinates of each DHS (150bp wide) observed in the genome in the
82 cell types, and field 4 holds an 82-element, comma-delimited vector of
normalized and scaled signal intensities (one per cell type) observed at the DHS
displayed in fields 1-3.


## Obtaining a list of "promoter DHSs."

Correlations are computed between "promoter DHSs" and non-promoter DHSs ("distal
DHSs") within some radius of each promoter DHS, usually a 500kb radius around
each promoter DHS. Defining a promoter DHS for a transcription start site (TSS)
in a gene is an extremely inexact science. What I do for each TSS is to call the
DHS that overlaps it "the promoter DHS."  When a TSS is not overlapped by a DHS,
I assign a nearby DHS to be the promoter DHS, slightly favoring a DHS upstream
of the TSS over a DHS downstream from the TSS. When no DHS is found near enough
a TSS, or no sufficiently strong DHS is found near enough a TSS, I do not assign
a promoter DHS to that TSS, and that TSS therefore gets dropped from subsequent
correlation analyses.


My script for assigning promoter DHSs to TSSs allows the user to tune the
parameters used as thresholds and tiebreakers. My preferred default values are
given below.

You can of course define promoter DHSs however you like within reason.

To define promoter DHSs, you of course need to start with a list of TSSs. You
can use any list you like. The one I use is in this file: `EH_v2_TxStarts.bed6`

These are TSSs from Gencode v19 for which a member of our team found RNA-seq
support.

This is the command I submitted to our cluster, to define promoter DHSs:

    qsub -cwd -N GetProms -S /bin/bash getPromoterDHSs \
        EH_v2_TxStarts.bed6 masterDHSsAndTagCounts_82_hg19.bed4 10000 2500 20 100 promOutfile.bed13

You can also run it directly from the command line, if desired; this one doesn't
take very long to execute (a couple minutes, I think).

That script uses the BEDOPS program closest-features, which under certain
circumstances will return "NA" for the "feature" (here, a DHS) closest to an
input element (here, a TSS). Moreover, the ordering on each line (TSS, then DHS)
needs to be reversed. So after the script runs to completion, the first thing we
need to do is reverse the order and filter out any "NA" entries, by running a
command like this:

    awk 'BEGIN{OFS="\t"}{if($13!="NA"){print $5,$6,$7,$4,".",$8,$1,$2,$3;}}' \
        promOutfile.bed13 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9

(The middle commands are necessary to remove any duplicate rows arising from
distinct lines in the TSS file whose first 4 columns are identical. If you use a
TSS file that doesn't contain any such rows, the middle commands would be
superfluous. I turn tab (`\t`) characters into '@' characters because I know
there are no `@` characters in the gene names, chromosome names, or
coordinates.)

Often, multiple positions within a few bp of one another get called as TSSs for
the same gene. Such TSSs will map to the same DHS. So next, we need to create a
file containing only the unique DHSs, by doing this:

    cut -f1-4 promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9 \
        | tr '\t' '@' \
        | sort -k 1b,1 \
        | uniq \
        | tr '@' '\t' \
        | sort-bed - \
        > promDHSsWithGeneNames_EHv2_82celltypes_hg19.bed4

Lastly, we create a file in which each row contains the promoter DHS, gene name,
and vector of normalized counts from the master list file, i.e., the vector that
will be used in the correlation calculations:

    bedmap --exact --delim "\t" --echo --echo-map-id \
        promDHSsWithGeneNames_EHv2_82celltypes_hg19.bed4 \
        masterDHSsAndTagCounts_82_hg19.bed4 \
        > promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5

## Obtaining a list of "distal DHSs."

"Distal DHSs" are all DHSs in the master list of DHSs that aren't promoter DHSs,
except:

* They should be limited to DHSs within a desired radius of promoters, e.g.,
  500kb
* DHSs very close to promoter DHSs should be excluded

So, do the following:

    cut -f1-3 promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5 \
        | uniq \
        | bedops -n -1 masterDHSsAndTagCounts_82_hg19.bed4 - \
        > nonprom_temp.bed4

This gives the DHSs which aren't promoter DHSs. Exclude DHSs that are, say, 1kb
from a promoter DHS:

    closest-features --shortest --dist nonprom_temp.bed4 \
        promDHSsAndTheirTSSs_EHv2_82celltypes_hg19.bed9 \
        | awk -F "|" '{dist=$3;if(dist<0){dist=-dist;}if(dist>=1000){print $1;}}' \
            > nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS_82celltypes.bed4

Finally, filter these DHSs to include only those within your target distance
from the nearest promoter DHS. Use the following command. The file it will
produce will not be easily readable like the file of promoter DHSs, but it will
be in the format expected by the program I use to compute correlations.
`--range 500000` restricts to distal DHSs within 500kb of promoter DHSs; you can
of course use a different radius if desired.

    bedmap --range 500000 --skip-unmapped --echo --echo-map \
    promDHSsWithGeneNamesAndTagCounts_EHv2_82celltypes_hg19.bed5 \
    nonPromoterDHSsAndTheirTagCounts_eachAtLeast1kbFromPromDHS_82celltypes.bed4 \
      > inputForCorrelationCalculations_500kb_82celltypes.bed


## Calculating correlations and storing them in files.

First, make sure you've followed the instructions at the outset in the first
section, "Preliminary Steps," to make an executable file on your system for the
correlation program.

Then, pick a threshold you're particularly interested in (e.g., 0.7), and supply
this and your preferred names for the output files to the following script, like
this:

    qsub -cwd -N GetCorrs500kb0.7 -S /bin/bash calcCorrelations.sh \
        inputForCorrelationCalculations_500kb_82celltypes.bed 0.7 corrs_promsFirst_EHv2_all_82celltypes_500kb.bed8 \
        corrs_promsFirst_EHv2_above0.7_82celltypes_500kb.bed8 \
        corrs_distalsFirst_EHv2_all_82celltypes_500kb.bed8 \
        corrs_distalsFirst_EHv2_above0.7_82celltypes_500kb.bed8

It will create 4 output files. One will contain all distal/promoter pairs and
their correlations (Pearson's _r_), in the format promoter, gene name, distal
DHS, _r_. Another will contain the subset of this, restricted to entries with
_r_ > your threshold. The other two files will be identical to these two, except
they'll be in the format distal DHS, gene name, promoter DHS, _r_. This might
give you 2 files you don't want, or 3 files you don't want, or you might find
uses for all 4. You can delete what you don't want, and use what you want.

Once these files are created, it's a good idea to delete the input file (here,
`inputForCorrelationCalculations_500kb_82celltypes.bed`), becuase it's really
large, no longer needed, and can be recreated again if desired.

That's all!
