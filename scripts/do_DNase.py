#!/usr/bin/env python

import sys
from multiprocessing import Pool
import subprocess

d_dir = sys.argv[1] # Directory where the fold change signal bigwigs for DNase signal are stored
o_dir = sys.argv[2] # Output directory

filenames = subprocess.check_output(["ls", d_dir]).split()

def wrapper(fname):
	#fname = DNASE.HepG2.fc.signal.bigwig
	cell = fname.split(".")[1]

	subprocess.call(
		" ".join([
			"bigWigAverageOverBed", 
			"-minMax", 
			d_dir+"/"+fname,
			"average_regions.bed",
			o_dir + "/" + cell + "_DNase.tsv"
		]), 
		shell = True
	)

p = Pool(len(filenames))
p.map(wrapper, filenames)
p.close()
p.join()
