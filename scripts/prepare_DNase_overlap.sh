#!/usr/bin/env bash

source activate DREAM-TFBS

data_dir='/cache/DREAM-tfbs-pred/'
out_dir=$data_dir'/out/'

ls /cache/DREAM-tfbs-pred/DNase/conservative/ | cut -f2 -d '.' | while read cell
  do
    bedtools intersect \
      -loj -a $data_dir'/annotations/regions/train_regions.blacklistfiltered.bed'  \
      -b  $data_dir'/DNase/conservative/DNASE.'$cell'.conservative.narrowPeak' \
      > $out_dir$cell'_DNASE_intersect_train.out' &
    
    bedtools intersect \
      -loj -a $data_dir'/annotations/regions/test_regions.blacklistfiltered.bed'  \
      -b  $data_dir'/DNase/conservative/DNASE.'$cell'.conservative.narrowPeak' \
      > $out_dir$cell'_DNASE_intersect_test.out' &
    
    bedtools intersect \
      -loj -a $data_dir'/annotations/regions/ladder_regions.blacklistfiltered.bed'  \
      -b  $data_dir'/DNase/conservative/DNASE.'$cell'.conservative.narrowPeak' \
      > $out_dir$cell'_DNASE_intersect_ladder.out' &
  done