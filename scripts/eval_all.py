#!/usr/bin/env python

import sys
import subprocess

data_dir = sys.argv[1]
out_dir = sys.argv[2]

f = open("TF_cellTypes.txt")
f.readline()

TFs = {}
all_cells = []
for line in f:
	fields = line.strip("\n").split("\t")
	TFs[fields[0]] = {
		"train":  fields[1].split(", "), 
		"L": fields[2].split(", "),  # ladder round
		"F":   fields[3].split(", ") # final test round
		}
	for Round in TFs[fields[0]]:
		for cell in Round:
			if cell not in all_cells: all_cells.append(cell)

command = """TF=%s; cell=%s; data_dir=%s; out_dir=%s;
if [ ! -e $out_dir'/%s.'$TF'.'$cell'.tab.gz' ]
  then
    if [ ! -e $out_dir'/%s.'$TF'.'$cell'.tab' ]
      then 
        python scripts/train_test_classifier.py \
        $data_dir'/annotations/regions/%s' \
        $TF \
        $out_dir'/'$TF'_PWM.tsv' \
        $out_dir'/'$TF'_StruM.tsv' \
        $out_dir'/'$cell'_DNase.tsv' \
        $out_dir'/RNA_vals.tsv' \
        $cell \
        $out_dir'/kmers.tsv' \
        $out_dir \
        $out_dir'/'$TF'_'*'_train.h5' \
        > $out_dir'/%s.'$TF'.'$cell'.tab'
  fi;
  gzip $out_dir'/%s.'$TF'.'$cell'.tab'
fi
"""

for round in ['L', 'F']:
	for TF in sorted(TFs.keys()):
		if sum( [1 for cell in TFs[TF]["F"] if cell != ""] ) == 0:
			continue
		if round == "L": region_file = "ladder_regions.blacklistfiltered.bed"
		else: region_file = "test_regions.blacklistfiltered.bed"
		for cell in TFs[TF][round]:
			if cell != "":
				subprocess.call(command % (TF, cell, data_dir, out_dir, round, round, region_file, round, round),
					shell=True)













