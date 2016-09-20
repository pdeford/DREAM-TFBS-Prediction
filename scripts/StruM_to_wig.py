#!/usr/bin/env python

"""
Usage: ./script.py TF TF_StruM.p TF_PWM.p shape_dir > TF_StruM.wig
"""

import sys
import numpy as np
import pickle
import subprocess
from multiprocessing import Pool
import bx.bbi.bigwig_file

TF = sys.argv[1]
StruM = pickle.load(open(sys.argv[2]))
PWM = pickle.load(open(sys.argv[3]))
shape_dir = sys.argv[4]
p = len(subprocess.check_output(["ls",shape_dir]).split())
k = len(StruM[0])/p

def norm_pdf(x,mu,var):
	l = (1./np.sqrt(2*np.pi*var))
	r = np.exp(-1*(x-mu)**2/(2*var))
	result = l*r
	#result = (1./np.sqrt(2*np.pi*var))*np.exp(-1*(x-mu)**2/(2*var))
	result += 10**-300
	return result

def score_strum(StruM,kmer_struc):
	return np.sum(np.log2(norm_pdf(kmer_struc,StruM[0],StruM[1]**2)))

def matchStrum(TF,StruM,kmer):
	kmer2 = rev_comp(kmer)
	p1 = score_strum(StruM,kmer)
	p2 = score_strum(StruM,kmer2)
	return max(p1,p2)

def rev_comp(struc):
	return np.hstack([struc[i:i+p] for i in range(0,len(struc),p)[::-1]])

def scoreStruM(TF,StruM,seq):
	return np.asarray( [matchStrum(TF,StruM,seq[i:i+(k)*p]) for i in range(0,len(seq)-(k+1)*p +1,p)] )

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

def get_sizes():
	sizes_path = "/cache/DREAM-tfbs-pred/annotations/hg19.chrom.sizes"
	sizes = {}
	for line in open(sizes_path):
		fields = line.split()
		sizes[fields[0]] = int(fields[1])
	return sizes

def lookup_seq_structure(shape_dir, chrom, start, end):
	data = []
	for f_name in subprocess.check_output(["ls","{}".format(shape_dir)]).split():
		bwh = bx.bbi.bigwig_file.BigWigFile(open(shape_dir + "/" + f_name, "rb"))
		row = bwh.get_as_array(chrom, start, end)
		row[np.isnan(row)] = 0.0
		data.append(row)

	data = np.vstack(data)
	return np.ravel(data)


def wrapper1(x):
	"""Load the structural sequence."""
	return lookup_seq_structure(shape_dir, chrom, x, x+100000)

def wrapper2(i):
	"""Score every position.."""
	return scoreStruM(TF,StruM,seq[i*p:(i+100000)*p])

sizes = get_sizes()
print "track type=wiggle_0"
for chrom_num in [str(i) for i in range(1,23)] + ["X"]:
	chrom = "chr" + chrom_num
	print >> sys.stderr, chrom
	print "fixedStep chrom={} start=1 step=1".format(chrom)
	stop = sizes[chrom]
	#seq = lookup_seq_structure(shape_dir, chrom, 1, stop)
	pool = Pool()
	seq = pool.map(wrapper1,range(0,stop,100000))
	pool.close()
	pool.join()
	seq = np.hstack(seq)[:stop*p]
	#scores = scoreStruM(TF,StruM,seq)
	pool = Pool()
	scores = pool.map(wrapper2, range(0,stop-k+1,100000-k+1))
	pool.close()
	pool.join()

	scores = np.hstack(scores)
	print "\n".join([str(x) for x in scores])






