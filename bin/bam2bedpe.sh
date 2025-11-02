#!/bin/bash

SDIR=$(dirname "$(readlink -f "$0")")

module load samtools
module load bedtools

BINS=data/raw/bins/hg38_20k_gz_enc_bins.bed
BTAG=$(basename ${BINS/_gz_*/})

BAM=$1

BASE=$(basename ${BAM/.bam/})
DIR=$(dirname $BAM)

ODIR=out/$DIR
mkdir -p $ODIR

samtools sort -n $BAM -@ 16 -m 1g \
    | samtools view -F 3084 -f 1 -b \
    | bedtools bamtobed -bedpe -i - \
    | tr ' ' '\t' \
    | gzip -9 - \
    > $ODIR/${BASE}.pe.bed.gz

Rscript $SDIR/src/map2bins.R $BINS $ODIR/${BASE}.pe.bed.gz

