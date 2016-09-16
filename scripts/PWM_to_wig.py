#!/usr/bin/env python

"""
Usage: ./script.py TF TF_PWM.p > TF_PWM.wig
"""

import sys
import numpy as np
import pickle
from multiprocessing import Pool

TF = sys.argv[1]
PWM = pickle.load(open(sys.argv[2]))

def rev_comp(seq):
	comp = dict(zip('ATGCN','TACGN'))
	return "".join([comp[n] for n in seq][::-1])

nuc_index = dict(zip("ACGTN",range(5)))
kmer_scores = {TF:{}}
def matchPWM(TF,PWM,kmer):
	k = PWM.shape[1]
	pwm = np.vstack([PWM,np.zeros([1,k])+0.25])
	if kmer in kmer_scores[TF]: 
		return kmer_scores[TF][kmer]
	else:
		kmer2 = rev_comp(kmer)
		p1 = np.product( [ pwm[nuc_index[n],j] for j,n in enumerate(kmer)] )
		p2 = np.product( [ pwm[nuc_index[n],j] for j,n in enumerate(kmer2)] )
		if k < 12:
			kmer_scores[TF][kmer] = kmer_scores[TF][kmer2] = max(p1,p2)
			return kmer_scores[TF][kmer]
		else:
			return max(p1,p2)

def scorePWM(TF,PWM,seq):
	k = PWM.shape[1]
	return np.asarray( [matchPWM(TF,PWM,seq[i:i+k]) for i in range(0,len(seq)-k +1)] )

def get_offsets():
	chrom_path = "/cache/DREAM-tfbs-pred/hg19.genome.fa"

	offsets = {}

	with open(chrom_path) as f:
		while 1:
			line = f.readline()
			pos = f.tell()
			if line.startswith(">"): offsets[line.strip()[1:]] = pos
			elif line.strip()  ==  "": break

	return offsets

def lookup_sequence(chrom, start=None, end=None, offsets=None):
	chrom_path = "/cache/DREAM-tfbs-pred/hg19.genome.fa"
	
	with open(chrom_path) as f:
		if offsets:
			pos = offsets[chrom]
		else:
			while 1:
				line = f.readline()
				if line.strip()[1:] == chrom:
					pos = f.tell()
					break
				elif line == "":
					return None
		if start is None and end is None:
			seq = ""
			f.seek(pos)
			while 1:
				line = f.readline().strip()
				if ">" in line: break
				seq += line
		elif not start is None and not end is None:
			f.seek( start + pos + (start//50) )
			lines = f.read(end - start + (end-start)//50 + 1)
			seq = "".join(lines.split("\n"))
		else: return None

	return seq

offsets = get_offsets()

def wrapper1(x):
	"""Scores."""
	return scorePWM(TF,PWM,seq[x:x+100000])

k = PWM.shape[1]
print "track type=wiggle_0"
for chrom_num in [str(i) for i in range(1,23)] + ["X"]:
	chrom = "chr" + chrom_num
	print "fixedStep chrom={} start=1 step=1".format(chrom)
	seq = lookup_sequence(chrom, offsets=offsets).upper()
	#scores = scorePWM(TF,PWM,seq)
	pool = Pool()
	scores = pool.map(wrapper1, range(0,len(seq)-k+1,100000-k+1))
	pool.close()
	pool.join()
	scores = np.hstack(scores)
	print "\n".join([str(x) for x in scores])






