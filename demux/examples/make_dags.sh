#!/bin/bash

set -euxo pipefail

IN_BCL="/mnt/nfs/bcldata/20260427_LH01056_0033_A23MTWJLT4"
SS="/mnt/groupdir/SUN-GI-GGSC-samplesheets/260427_samplesheet_PLUS25A.csv"
OUT_DIR="."

SNAKEMAKE="snakemake -c 4 -s ../../workflow/Snakefile --workflow-profile ../../../resources/profile/PROD.profile.yaml --forceall --config bcl=$IN_BCL sample_sheet=$SS out_dir=$OUT_DIR"

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
