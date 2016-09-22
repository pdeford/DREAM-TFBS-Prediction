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

command = """if [ ! -e $out_dir'/%s.'$TF'.'$cell'.tab' ]
  then
    TF=%s; cell=%s; data_dir=%s; out_dir=%s; python scripts/train_test_classifier.py \
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
"""

#for TF in TFs:
for TF in ["ARID3A"]:
	for round in ["L","F"]:
		if round == "L": region_file = "ladder_regions.blacklistfiltered.bed"
		else: region_file = "test_regions.blacklistfiltered.bed"
		for cell in TFs[TF][round]:
			if cell != "":
				subprocess.call(command % (round, TF, cell, data_dir, out_dir, region_file, round, round),
					shell=True)













