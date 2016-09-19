#!/usr/bin/env python

"""
Given input ChIP data for TFs, run `meme-chip` to learn 
PWM, and use `FIMO` output to learn StruM. Pickle the 
results as a np.array

Usage: cat ChIP.for.TF.*.bed | script.py TF out_directory shape_file_directory
"""

############################################################
# IMPORT NECESSARY MODULES                                 #
############################################################

from __future__ import division
import subprocess
import numpy as np
import pickle
import sys

#==========================================================#
# READ COMMAND LINE ARGUMENTS

TF = sys.argv[1]
chip = sys.stdin
outdir = sys.argv[2]
shapedir = sys.argv[3]

subprocess.call("mkdir %s/%s_meme" % (outdir,TF),shell=True)

############################################################
# DEFINE FUNCTIONS                                         #
############################################################

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

############################################################
# TAKE CHIP.BED AND CONVERT TO FASTA, HG19                 #
############################################################
print "Converting .bed to .fa"

offsets = get_offsets()

subprocess.call("mkdir %s" % outdir, shell=True)

TF_seq_fa = "%s/%s.fa" % (outdir, TF)
with open(TF_seq_fa, "wb") as g:
	for line in chip:
		fields = line.strip().split()
		
		chrom =     fields[0]
		start = int(fields[1])
		end   = int(fields[2])

		g.write( ">%s:%d-%d\n" % (chrom,start,end,) )
		seq = lookup_sequence(chrom,start,end,offsets)
		for i in range(0,len(seq),50): g.write(seq[i:i+50].upper() + "\n")


############################################################
# RUN MEME ON SEQUENCES TO GET PWM AND MATCHES             #
############################################################
print "Running MEME-ChIP on training sequences"

subprocess.call("mkdir output",shell=True)

subprocess.call(
	"meme-chip -oc %s/%s_meme -meme-nmotifs 1 -dreme-m 1  %s" % (
		outdir,
		TF, 
		TF_seq_fa
		),
	shell=True
	)

############################################################
# LOAD PWM MEME OUTPUT                                     #
############################################################
print "Reading PWM from MEME output"

with open("%s/%s_meme/meme_out/meme.txt" % (outdir,TF)) as f:
	nucs = "ACGT"
	PWM  = []
	start2 = False
	for line in f:
		if "Motif 1 position-specific probability matrix" in line:
			start2 = True
			count = -2
			continue
		if start2:
			count += 1
			if count > 0:
				if "--------" in line:
					start2 = False
					continue
				if line.strip() != "":
					A,C,G,T = [float(x) for x in line.split()]
					PWM.append([A,C,G,T])

	PWM = np.asarray(PWM).T + 0.01

############################################################
# LOAD STRUCTURE TRAINING SEQUENCES FROM FIMO OUTPUT       #
############################################################
print "Reading FIMO output"

with open("%s/%s_meme/fimo_out_1/fimo.txt" % (outdir,TF)) as f:
	f.readline()

	hits = []

	for line in f:
		fields = line.split()
		chrom  = fields[1]
		start  = int(fields[2])
		end    = int(fields[3])
		hits.append([chrom,start,end])

############################################################
# EXTRACT SHAPE AND DNAse SIGNAL FOR REGIONS               #
############################################################
print "Learning StruMs from FIMO output"



StruM = [[], []]

for f_name in subprocess.check_output(["ls","{}".format(shapedir)]).strip().split():
	if f_name == "": continue
	print "Pulling regions from {}".format(f_name)
	data = np.asarray(
		[
			[
				float(x) for x in subprocess.check_output(
						["bigWigSummary", "{}/{}".format(shapedir,f_name), chrom, "{}".format(start), "{}".format(end+1), "{}".format(end-start+1), "-type=mean"]
					).split()
			] for chrom,start,end in hits
		]
	)

	StruM[0].append(np.average(data,axis=0))
	StruM[1].append(np.std(data,axis=0))

StruM = [np.ravel(np.vstack(x)) for x in StruM]

############################################################
# PICKLE THE OUTPUT                                        #
############################################################

pickle.dump(PWM, open("{}/{}_PWM.p".format(outdir,TF),"wb"))
pickle.dump(StruM, open("{}/{}_StruM.p".format(outdir,TF),"wb"))

