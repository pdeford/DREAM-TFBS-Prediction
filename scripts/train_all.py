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
		"ladder": fields[2].split(", "),
		"test":   fields[3].split(", ")
		}
	for Round in TFs[fields[0]]:
		for cell in Round:
			if cell not in all_cells: all_cells.append(cell)


#for TF in TFs:
for TF in ["ATF7"]:
	for cell in TFs[TF]["train"]:
		subprocess.call(

		"""TF=%s; cell=%s; data_dir=%s; out_dir=%s; python scripts/build_train_array.py \
		$data_dir'/annotations/labels/'$TF'.train.labels.tsv' \
		$TF \
		$out_dir'/'$TF'_PWM.tsv' \
		$out_dir'/'$TF'_StruM.tsv' \
		$out_dir'/'$cell'_DNase.tsv' \
		$out_dir'/RNA_vals.tsv' \
		$cell \
		$out_dir'/kmers.tsv' \
		$out_dir""" % (TF, cell, data_dir, out_dir),
			shell=True)


