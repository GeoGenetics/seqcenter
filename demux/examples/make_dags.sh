#!/bin/bash

set -euxo pipefail

IN_BCL="/mnt/nfs/bcldata/20260810_LH01056_0075_A23JC7CLT3/"
SS="/mnt/groupdir/SUN-GI-GGSC-samplesheets/260810_samplesheet_PLUS10A.csv"
OUT_DIR=`pwd`

SNAKEMAKE="snakemake --snakefile ../../workflow/Snakefile --workflow-profile ../../resources/profile/PROD.profile.yaml --forceall --config bcl=$IN_BCL sample_sheet=$SS out_dir=$OUT_DIR idx_list=resources/IndexInformation_20260819JBT.tsv"

cd bcl_convert
$SNAKEMAKE demux=yes --dryrun
$SNAKEMAKE demux=yes --rulegraph | dot -Tsvg > rulegraph.svg
$SNAKEMAKE demux=yes --filegraph | dot -Tsvg > filegraph.svg
$SNAKEMAKE demux=yes --dag | dot -Tsvg > dag.svg
cd ../


cd copy
$SNAKEMAKE demux=no --dryrun
$SNAKEMAKE demux=no --rulegraph | dot -Tsvg > rulegraph.svg
$SNAKEMAKE demux=no --filegraph | dot -Tsvg > filegraph.svg
$SNAKEMAKE demux=no --dag | dot -Tsvg > dag.svg
cd ../
