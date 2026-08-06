#!/bin/bash

shopt -s extglob
set -euo pipefail

module load pixi/0.65.0

IN_BCL=$1; shift
SS=$1; shift
OUT_DIR=$1; shift

# analysis_n=[last]
# resume=False (True/False)
# suffix=""
pixi run --manifest-path /projects/ggsc/apps/seqcenter/demux/ snakemake -c 4 -s /projects/ggsc/apps/seqcenter/demux/workflow/Snakefile --workflow-profile /projects/ggsc/apps/resources/profile/PROD.profile.yaml --config bcl=$IN_BCL sample_sheet=$SS out_dir=$OUT_DIR $@
