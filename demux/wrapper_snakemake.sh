#!/bin/bash

set -euo pipefail
export XDG_CACHE_HOME=/projects/ggsc/scratch
export PIXI_CACHE_DIR=/tmp/pixi/cache

module load pixi/0.77.1

BASEDIR=$(dirname `realpath $0`)
IN_BCL=`realpath --canonicalize-existing --no-symlinks $1`; shift
SS=`realpath --canonicalize-existing --no-symlinks $1`; shift
OUT_DIR=`realpath --canonicalize-existing --no-symlinks $1`; shift

# Define workdir (for temp/, logs/, and .snakemake/)
WORKDIR=/projects/ggsc/scratch/demux/$USER
mkdir -p $WORKDIR

# Set conda channel priority
pixi run --as-is --manifest-path $BASEDIR conda config --set channel_priority strict

# Run workflow
env --chdir=$WORKDIR pixi run --as-is --manifest-path $BASEDIR snakemake --snakefile $BASEDIR/workflow/Snakefile --workflow-profile /projects/ggsc/apps/seqcenter/demux/resources/profile.yaml --config bcl=$IN_BCL sample_sheet=$SS out_dir=$OUT_DIR $@
