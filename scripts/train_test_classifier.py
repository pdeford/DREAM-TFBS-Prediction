#!/usr/bin/env python

from __future__ import division
import sys
import numpy as np
import pickle

regions = open(sys.argv[1]) # data/annotations/labels/TF.train.labels.tsv
TF = sys.argv[2] # Duh
pwm_file = open(sys.argv[3]) # output/TF_PWM.tsv
strum_file = open(sys.argv[4]) # output/TF_StruM.tsv
dnase_file = open(sys.argv[5]) # output/DNASE_cell.tsv
rna_file = open(sys.argv[6]) # output/RNA_vals.tsv
train_cell = sys.argv[7] # Duh
kmer_file = open(sys.argv[8]) # output/kmers.tsv
training_arrays = np.vstack([pickle.load(open(x)) for x in sys.argv[9:]])

#============================================================
# PREPROCESS ARRAY
print >> sys.stderr, "Loading training data"

Y = training_arrays[:,0]
training_arrays = np.delete(training_arrays,0,1)
X = training_arrays

clean_avg = np.average(X, axis=0)
clean_std = np.std(X,axis=0)

X = (X-clean_avg)/clean_std

#============================================================
# TRAIN THE MODEL
print >> sys.stderr, "Picking model"

from sklearn.svm import SVC
from sklearn.linear_model import LogisticRegression as logit
from sklearn.ensemble import RandomForestClassifier as RFC
from sklearn.naive_bayes import GaussianNB as GNB

from sklearn.metrics import roc_auc_score, average_precision_score
from sklearn.cross_validation import StratifiedKFold


kfold = 5
skf = StratifiedKFold(Y, n_folds=kfold, shuffle=True, random_state=1104)

# Learn params
svc_clf = GridSearchCV(
	SVC(probability=True, n_jobs=-1),
	param_grid={
		'kernel':['rbf','poly'], 
		'C':[0.01, 0.1, 1.0, 10, 100],
		'degree': [3]
		},
	cv=skf,
	)
svc_clf.fit(X,Y)
svc_clf = SVC().set_params(**svc_clf.get_params(deep=True))

log_clf = GridSearchCV(
	logit(solver='sag', n_jobs=-1),
	param_grid={
		'penalty':['l1','l2'],
		'C':[0.01, 0.1, 1.0, 10, 100],
		},
	cv=skf,
	)
log_clf.fit(X,Y)
log_clf = logit().set_params(**log_clf.get_params(deep=True))

rfc_clf = RFC(n_estimators=20)	
gnb_clf = GNB()

# Cross validation to pick the appropriate classifier
clfs = [svc_clf, log_clf, rfc_clf, gnb_clf]
scores = []
for clf in clfs:
	score = 0
	for train_index, test_index in skf:
		clf.fit(X[train_index],Y[train_index])
		score += average_precision_score(Y[test_index],clf.predict_proba(X[test_index]))
	scores.append(score)

n = np.argmax(scores)
best_clf = clfs[n]
best_clf.fit(X)
print >> sys.stderr, "Best classifier:\n\tarPRC: {}\n\tClassifier: {}\n\tParameters: {}".format(scores[n]/kfold, best_clf, best_clf.get_params())

del X, Y, training_arrays

#============================================================
# BUILD TEST ARRAY

print >> sys.stderr, "Loading test data"

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


header = rna_file.readline().strip().split()
cell_columns = [i for i,c in enumerate(header) if c==train_cell]

data = []

chroms = []
for line in regions:

	fields = line.strip().split()
	chrom, start, end = fields[0], int(fields[1]), int(fields[2])
	mid = (start+end)/2

		for rna_line in rna_file:
		rna_fields = rna_line.strip().split()
		if rna_fields[0] == chrom:
			if (int(rna_fields[1])+int(rna_fields[2]))/2 == mid:
				#RNA = [float(x) for x in rna_fields[3:]]
				RNA = [float(rna_fields[i]) for i in cell_columns]
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

#============================================================
# PROCESS THE DATA
print >> sys.stderr, "Scaling data"


data = np.asarray(data)
data = (data-clean_avg)/clean_std

#============================================================
# SCORE EACH POSITION
print >> sys.stderr, "Scoring"

Y = best_clf.predict_proba(data)
del data

regions.seek(0)
count = 0
for line in regions:
	print line.strip() + "\t{}".format(Y[count])
	count += 1

