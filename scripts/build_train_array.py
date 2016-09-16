#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
import pickle

"""
Usage:
python build_train_array.py \
  $data_dir'/annotations/labels/'$TF'.train.labels.tsv' \
  $TF \
  $out_dir'/'$TF'_PWM.tsv' \
  $out_dir'/'$TF'_StruM.tsv' \
  $out_dir'/'$cell'_DNase.tsv' \
  $out_dir'/RNA_vals.tsv' \
  $cell \
  $out_dir'/kmers.tsv'
"""


training_chip = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
TF = sys.argv[2] # Duh
pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
train_cell = sys.argv[7] # Duh
kmer_file = open(sys.argv[8]) # output/kmers.tsv


def calc_foot(scores):
	return [(DNASE[i-2]+DNASE[i+2])/2-DNASE[i] * scores[i] for i in range(2,len(DNASE)-2)]

def update_data(file_object, column, object_to_update, counter, column_end=None):
	if column_end is None: slicer = lambda x:float(x[column])
	elif column_end is True: slicer = lambda x:[float(f) for f in  x[column:]]
	else: slicer = lambda x:[float(f) for f in x[column:column_end]]

	count = 0
	for line in file_object:
		count += 1
		fields = line.strip().split()
		object_to_update.pop(0)
		object_to_update.append(slicer(fields))
		if counter == counter: break

def populate(file_object, column, interval, size, column_end=None):
	if column_end is None: slicer = lambda x:float(x[column])
	elif column_end is True: slicer = lambda x:[float(f) for f in x[column:])
	else: slicer = lambda x:[float(f) for f in x[column:column_end]]

	l = []
	for line in file_object:
		fields = line.strip().split()
		position = fields[0].split(":")
		if position[0] == chrom:
			if int(position.split("-")) == start + len(l)*size-(interval-200)/2:
				l.append(slicer(fields))
				if len(l) == interval/size: break
	return l



header = training_chip.readline().strip().split()
for i in range(len(header)):
	if header[i] == train_cell:
		break
cell_column = i


data = []

chroms = []
for line in training_chip:

	fields = line.strip().split()
	chrom, start, end, bound = fields[0], int(fields[1]), int(fields[2]), fields[cell_column] 
	mid = (start+end)/2
	
	if bound == "A": continue
	if bound == "B": Y = 1
	else: Y = 0


	for rna_line in rna_file:
		rna_fields = rna_line.strip().split()
		if rna_fields[0] == chrom:
			if (int(rna_fields[1])+int(rna_fields[2]))/2 == mid:
				RNA = [float(x) for x in rna_fields[3:]]
				break

	if chrom not in chroms:
		chroms.append(chrom)
		#for f in [pwm_file, strum_file, dnase_fil, rna_file]: f.seek(0)
	
		DNASE  = populate(dnase_file, 4, 500, 25)
		PWM    = populate(pwm_file, 7, 500, 25)
		STRUM  = populate(strum_file, 7, 500 25)
		KMERS  = populate(kmer_file, 1, 500, 25,True)
		P_FOOT = calc_foot(PWM)
		S_FOOT = calc_foot(STRUM)
	else:
		interval = start-last_start
		mult25 = interval/25
		
		update_data(dnase_file, 4, DNASE, mult25)
		update_data(pwm_file, 7, PWM, mult25) 
		update_data(strum_file, 7, STRUM, mult25)
		update_data(kmer_file, 1, KMERS, mult25, True)
		
		P_FOOT = calc_foot(PWM)
		S_FOOT = calc_foot(STRUM)

	row = Y + RNA + DNASE + PWM + STRUM + P_FOOT + S_FOOT+ list(np.sum(KMERS, axis=0))
	data.append(row)
	
	last_start = start

pickle.dump(data,open("{}_{}_train.p".format(TF,train_cell),"wb"))


