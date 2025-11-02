#!/bin/bash

rm -vrf out
./bin/bam2bedpe.sh data/raw/test/bams/ALT3891/ALT3891_A_1056/b188/ALT3891_A_1056_b188___MD.bam
zcat out/20k/data/raw/test/bams/ALT3891/ALT3891_A_1056/b188/ALT3891_A_1056_b188___MD.pe.map.gz  | md5sum -
zcat data/raw/pe.map/20k/1Cell/ALT3891/ALT3891_A_1056/b188/ALT3891_A_1056_b188___MD.pe.map.gz | md5sum -

if [ "$FULL" == "Yes" ]; then
  ./bin/bam2bedpe.sh data/raw/test/bams/ALT3891/ALT3891_A_1048/b63/ALT3891_A_1048_b63___MD.bam
  zcat data/raw/pe.map/20k/1Cell/ALT3891/ALT3891_A_1048/b63/ALT3891_A_1048_b63___MD.pe.map.gz | md5sum -
  zcat out/20k/data/raw/test/bams/ALT3891/ALT3891_A_1048/b63/ALT3891_A_1048_b63___MD.pe.map.gz | md5sum -
fi
