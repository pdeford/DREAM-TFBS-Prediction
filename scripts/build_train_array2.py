#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
import cPickle as pickle
from multiprocessing import Pool

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
  $out_dir'/kmers.tsv' \
  $out_dir
"""

print >> sys.stderr, "Reading arguments"

training_chip = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
TF = sys.argv[2] # Duh
pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
train_cell = sys.argv[7] # Duh
kmer_file = open(sys.argv[8]) # output/kmers.tsv
out_dir = sys.argv[9]


def calc_foot(scores):
	return [(DNASE[i-2]+DNASE[i+2])/2-DNASE[i] * scores[i] for i in range(2,len(DNASE)-2)]



print >> sys.stderr, "Ready, set..."

header = training_chip.readline().strip().split()
for i in range(len(header)):
	if header[i] == train_cell:
		break
cell_column = i

header = rna_file.readline().strip().split()
cell_columns = [i for i,c in enumerate(header) if c==train_cell]

data = []

print >> sys.stderr, "Go!"

def wrapper(line):
	def populate(file_object, column, interval, size, column_end=None):
		if column_end is None: slicer = lambda x:float(x[column])
		elif column_end is True: slicer = lambda x:[float(f) for f in x[column:]]
		else: slicer = lambda x:[float(f) for f in x[column:column_end]]

		l = []
		for line in file_object:
			fields = line.strip().split()
			position = fields[0].split(":")
			if position[0] == chrom:
				if int(position[1].split("-")[0]) == start + len(l)*size-(interval-200)/2:
					l.append(slicer(fields))
					if len(l) == interval/size: break
		return l

	pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
	strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
	dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
	rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
	kmer_file = open(sys.argv[8]) # output/kmers.tsv

	fields = line.strip().split()
	chrom, start, end, bound = fields[0], int(fields[1]), int(fields[2]), fields[cell_column] 
	mid = (start+end)/2
	
	if bound == "A": return None
	if bound == "B": Y = 1
	else: Y = 0

	for rna_line in rna_file:
		rna_fields = rna_line.strip().split()
		if rna_fields[0] == chrom:
			if int(rna_fields[1]) == mid:
				#RNA = [float(x) for x in rna_fields[3:]]
				RNA = [float(rna_fields[i]) for i in cell_columns]
				break

	DNASE  = populate(dnase_file, 4, 500, 25)
	PWM    = populate(pwm_file, 7, 500, 25)
	STRUM  = populate(strum_file, 7, 500, 25)
	KMERS  = populate(kmer_file, 1, 500, 25,True)
	P_FOOT = calc_foot(PWM)
	S_FOOT = calc_foot(STRUM)

	print "From here"
	row = [Y] + RNA + DNASE + PWM + STRUM + P_FOOT + S_FOOT+ list(np.sum(KMERS, axis=0))
	print "to HERE!"

	return row


while True:
	breaker = False
	lines = [training_chip.readline() for i in range(100000)]
	poppers = [i for i,l in enumerate(lines) if l == ""][::-1]
	for i in poppers: 
		lines.pop(i)
		breaker = True
	pool = Pool()
	sub_data = pool.map(wrapper, lines)
	pool.close()
	pool.join()

	for row in sub_data: 
		if row is None:
			pass
		else:
			data.append(row)

	if breaker: break


pickle.dump(data,open(out_dir+"/{}_{}_train.p".format(TF,train_cell),"wb"))


