#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
from multiprocessing import Pool, cpu_count
import subprocess
import tables
import glob

"""
Usage:

out_dir=/cache/DREAM-tfbs-pred/out/
data_dir=/cache/DREAM-tfbs-pred/
TF=ATF7
cell=K562

python build_train_array2.py \
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

pwm_file_names = glob.glob(out_dir + "/*_PWM.tsv")
strum_file_names = glob.glob(out_dir + "/*_StruM.tsv")

wc_l = sum([1 for line in training_chip])
training_chip.seek(0)

print >> sys.stderr, "Ready, set..."

header = training_chip.readline().strip().split()
for i in range(len(header)):
	if header[i] == train_cell:
		break
cell_column = i

positions = [training_chip.tell()]
n_cores = cpu_count()
chunk_size = wc_l//n_cores + 1
for i in range(n_cores):
	for i in range(chunk_size): training_chip.readline()
	positions.append(training_chip.tell())


header = rna_file.readline().strip().split()
cell_columns = [i for i,c in enumerate(header) if c==train_cell]

def wrapper(start_stop):

	f_start = start_stop[0]
	f_stop = start_stop[1]


	training_chip = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
	pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
	pwm_files = [open(x) for x in pwm_file_names]
	strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
	strum_files = [open(x) for x in strum_file_names]
	dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
	rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
	kmer_file = open(sys.argv[8]) # output/kmers.tsv
	
	training_chip.seek(f_start)
	f_pos = f_start

	data = []
	Y = []

	chroms = []

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

	while f_pos < f_stop:
		line = training_chip.readline()
		f_pos = training_chip.tell()

		fields = line.strip().split()
		chrom, start, end, bound = fields[0], int(fields[1]), int(fields[2]), fields[cell_column] 
		mid = (start+end)/2
		
		if bound == "A": continue
		if bound == "B": y = 1
		else: y = 0


		for rna_line in rna_file:
			rna_fields = rna_line.strip().split()
			if rna_fields[0] == chrom:
				if int(rna_fields[1]) == mid:
					#RNA = [float(x) for x in rna_fields[3:]]
					RNA = [float(rna_fields[i]) for i in cell_columns]
					break

		if chrom not in chroms:
			#print >> sys.stderr, "..." + chrom
			chroms.append(chrom)
			for f in strum_files + pwm_files + [pwm_file, strum_file, dnase_file, rna_file]: f.seek(0)
		
			DNASE  = populate(dnase_file, 4, 500, 25)
			PWM    = populate(pwm_file, 7, 500, 25)
			PWM2   = [populate(p_file, 7, 500, 25) for p_file in pwm_files]
			STRUM  = populate(strum_file, 7, 500, 25)
			STRUM2   = [populate(s_file, 7, 500, 25) for s_file in strum_files]
			KMERS  = populate(kmer_file, 1, 500, 25,True)
			P_FOOT = calc_foot(PWM)
			S_FOOT = calc_foot(STRUM)
		else:
			interval = start-last_start
			mult25 = interval/25
			
			update_data(dnase_file, 4, DNASE, mult25)
			update_data(pwm_file, 7, PWM, mult25)
			for i,pwm_mini in enumerate(PWM2):
				update_data(pwm_files[i], 7, pwm_mini, mult25)
			update_data(strum_file, 7, STRUM, mult25)
			for i,strum_mini in enumerate(STRUM2):
				update_data(strum_files[i], 7, strum_mini, mult25)
			update_data(kmer_file, 1, KMERS, mult25, True)
			
			P_FOOT = calc_foot(PWM)
			S_FOOT = calc_foot(STRUM)

		row = RNA + DNASE + PWM + STRUM + P_FOOT + S_FOOT+ list(np.sum(KMERS, axis=0)) + \
			[np.max(DNASE), np.max(PWM), np.max(STRUM), np.max(P_FOOT), np.max(S_FOOT)] + \
			[np.max(tf) for tf in PWM2] + [np.max(tf) for tf in STRUM2]

		Y.append(y)
		data.append(row)
		
		last_start = start


	data = np.array(data, dtype=np.float)
	Y = np.array(Y, dtype=np.int)
	name = out_dir+"/{}_{}_inter_{}.h5".format(TF,train_cell,f_start)
	h5file = tables.open_file(name, mode='w')
	#h5file.create_group(h5file.root, "data", "data")
	h5file.create_array(h5file.root, 'X', data, 'data')
	h5file.create_array(h5file.root, 'Y', Y, 'data')
	h5file.close()
	return name

print >> sys.stderr, "Go!"

pool = Pool(n_cores)
data_names = pool.map(wrapper, [(positions[i],positions[i+1]) for i in range(len(positions)-1)])
pool.close()
pool.join()

print >> sys.stderr, "Done, compiling"
h5file = tables.open_file(out_dir+"/{}_{}_train.h5".format(TF,train_cell), mode = 'w')
data_group = h5file.create_group(h5file.root, "data", "{} data".format(TF))
for i,name in enumerate(data_names):
	h5file2 = tables.open_file(name, mode="r")
	x = h5file2.root.X[:]
	y = h5file2.root.Y[:]
	if i == 0:
		X = h5file.create_earray(data_group, 'X', tables.Float32Atom(), [0, x.shape[1]], "matrix values", expectedrows=52000000)
		Y = h5file.create_earray(data_group, 'Y', tables.IntAtom(), [0,], "bound values", expectedrows=52000000)
	X.append(x)
	Y.append(y)
	h5file2.close()
	subprocess.call("rm {}".format(name), shell=True)
	

print >> sys.stderr, "Done. Closing"
h5file.close()

