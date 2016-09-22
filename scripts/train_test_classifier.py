#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
import pandas as pd
import cPickle as pickle
from multiprocessing import Pool, cpu_count
import subprocess
import random
import tables
import os.path

regions = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
TF = sys.argv[2] # Duh
pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
test_cell = sys.argv[7] # Duh
kmer_file = open(sys.argv[8]) # output/kmers.tsv
out_dir = sys.argv[9]
training_arrays = sys.argv[10:]

n_cores = cpu_count()

if os.path.isfile(out_dir + "/%s_clf.p" % TF):
	best_clf = pickle.load(open(out_dir + "/%s_clf.p" % TF,"rb"))
else:
	#============================================================
	# PREPROCESS ARRAY
	print >> sys.stderr, "Loading training data"

	N1 = 10000
	N2 = 1000000 # 1Mil positive examples (a full cell type's data's worth) and 10Mil neg (1/5 total)
	n1 = int(N1/len(training_arrays))
	n2 = int(N2/len(training_arrays))

	X_train = []
	Y_train = []
	X_test = []
	Y_test = []

	for name in training_arrays:
		f = tables.open_file(name, mode='r')
		y = f.root.data.Y[:]
		x = f.root.data.X

		all_indices = np.asarray(range(y.shape[0]))
		pos_sample = all_indices[y==1]
		neg_sample = all_indices[y==0]

		np.random.shuffle(pos_sample)
		np.random.shuffle(neg_sample)

		training_indices = list(np.hstack([pos_sample[:n1], neg_sample[:10*n1]]))
		full_train = list(np.hstack([pos_sample[n1:n1+n2], neg_sample[10*n1:10*(n1+n2)]]))

		X_test.append(x[training_indices,:])
		Y_test.append(y[training_indices])
		X_train.append(x[full_train,:])
		Y_train.append(y[full_train])
		f.close()
		del y,x,all_indices,pos_sample,neg_sample,training_indices,full_train

	X_train = np.vstack(X_train)
	X_test = np.vstack(X_test)
	Y_train = np.hstack(Y_train)
	Y_test = np.hstack(Y_test)

	print >> sys.stderr, "Normalizing data"
	clean_avg = np.average(X_train, axis=0)
	clean_std = np.std(X_train,axis=0)

	def normalize(array):
		return (array-clean_avg)/clean_std

	X_test = normalize(X_test)
	X_train = normalize(X_train)

	#============================================================
	# TRAIN THE MODEL
	print >> sys.stderr, "Training model"

	from sklearn.ensemble import RandomForestClassifier as RFC
	from sklearn.metrics import roc_auc_score, average_precision_score

	best_clf = RFC(n_jobs=-1,n_estimators=50)
	best_clf.fit(X_train, Y_train)
	Y2 = best_clf.predict_proba(X_test)[:,1]
	print >> sys.stderr, best_clf
	print >> sys.stderr, "auPRC:", average_precision_score(Y_test, Y2)
	print >> sys.stderr, "auROC:", roc_auc_score(Y_test, Y2)

	pickle.dump(best_clf, open(out_dir + "/%s_clf.p" % TF,"wb"))

	training_arrays.close()
	del X_train, X_test, Y_train, Y_test

#============================================================
# BUILD TEST ARRAY

print >> sys.stderr, "Loading test data"

wc_l = sum([1 for line in regions])
regions.seek(0)
positions = [regions.tell()]
n_cores = cpu_count()
chunk_size = wc_l//n_cores + 1
for i in range(n_cores):
	for i in range(chunk_size): regions.readline()
	positions.append(regions.tell())


def wrapper(start_stop):

	regions = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
	pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
	strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
	dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
	rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
	kmer_file = open(sys.argv[8]) # output/kmers.tsv

	f_start = start_stop[0]
	f_stop = start_stop[1]
	f_pos = f_start
	regions.seek(f_start)

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


	header = rna_file.readline().strip().split()
	cell_columns = [i for i,c in enumerate(header) if c==test_cell]

	data = []
	chroms = []

	while f_pos < f_stop:
		line = regions.readline()
		f_pos = regions.tell()
		fields = line.strip().split()
		chrom, start, end = fields[0], int(fields[1]), int(fields[2])
		mid = (start+end)/2

		for rna_line in rna_file:
			rna_fields = rna_line.strip().split()
			if rna_fields[0] == chrom:
				if int(rna_fields[1]) == mid:
					#RNA = [float(x) for x in rna_fields[3:]]
					RNA = [float(rna_fields[i]) for i in cell_columns]
					break

		if chrom not in chroms:
			chroms.append(chrom)
			for f in [pwm_file, strum_file, dnase_file, rna_file]: f.seek(0)
		
			DNASE  = populate(dnase_file, 4, 500, 25)
			PWM    = populate(pwm_file, 7, 500, 25)
			STRUM  = populate(strum_file, 7, 500, 25)
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

		row = RNA + DNASE + PWM + STRUM + P_FOOT + S_FOOT+ list(np.sum(KMERS, axis=0))
		data.append(row)
		
		last_start = start

	name = out_dir+"/{}_{}_inter_{}.h5".format(TF,test_cell,f_start)
	data = np.array(data, dtype=np.float)
	h5file = tables.open_file(name, mode='w')
	h5file.create_array(h5file.root, 'X', data, 'data')
	h5file.close()
	return name

pool = Pool(n_cores)
data_names = pool.map(wrapper, [(positions[i],positions[i+1]) for i in range(len(positions)-1)])
pool.close()
pool.join()

#============================================================
# PROCESS THE DATA
print >> sys.stderr, "Scaling data"

#data = np.vstack([pd.HDFStore(name)['data'] for name in data_names])
data = np.vstack([ tables.open_file(name).root.X[:] for name in data_names ])
data = np.asarray(data)
data = (data-clean_avg)/clean_std
for name in data_names: subprocess.call("rm %s" % name, shell=True)

#============================================================
# SCORE EACH POSITION
print >> sys.stderr, "Scoring"

Y = best_clf.predict_proba(data)[:,1]
del data

regions.seek(0)
count = 0
for line in regions:
	print line.strip() + "\t{}".format(Y[count])
	count += 1

