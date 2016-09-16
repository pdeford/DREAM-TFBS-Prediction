#!/usr/bin/env python

"""
Count the kmers in each sequence in the regions file. Count 
reverse strand as well.

Usage: ./script.py average_regions.bed k
"""

import sys
from multiprocessing import Pool

regions = open(sys.argv[1])
k = int(sys.argv[2])


def rev_comp(seq):
	comp = dict(zip('ATGCN','TACGN'))
	return "".join([comp[n] for n in seq][::-1])

def kmer_maker(k):
	nucs = "ATGC"
	kmers = [""] * len(nucs)**k
	for i in range(k):
		running = 0
		while running + 1 < len(kmers):
			for n in nucs:
				for j in range(len(nucs)**(k-(i+1))):
					kmers[j+running] += n
				running += len(nucs)**(k-(i+1))
	return kmers

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

offsets = get_offsets()

kmers = kmer_maker(k)
template = dict([(kmer,0) for kmer in kmers])

def wrapper(line):
	fields = line.split()
	chrom,start,stop = fields[0], int(fields[1]), int(fields[2])
        seq = lookup_sequence(chrom, start, stop+k, offsets).upper() 
	kmer_counts = template.copy()
	for i in range(len(seq)-k+1):
		kmer = seq[i:i+k]
		if "N" in kmer: continue
		kmer_counts[kmer] += 1
		kmer_counts[rev_comp(kmer)] += 1
	out = fields[3:4] + [str(kmer_counts[x]) for x in kmers]
	return "\t".join(out)

lines = []
while True:
	line = regions.readline().strip()
	if line == "": break
	lines.append(line)
	if len(lines) == 100000:
		pool = Pool()
		counts = pool.map(wrapper,lines)
		pool.close()
		pool.join()
		del lines
		lines = []
		for nl in counts: print nl
if len(lines) > 0:
	pool = Pool()
	counts = pool.map(wrapper,lines)
	pool.close()
	pool.join()
	del lines
	lines = []
	for nl in counts: print nl


#for line in regions:
#	fields = line.strip().split()
#	chrom,start,stop = fields[0], int(fields[1]), int(fields[2])
#	seq = lookup_sequence(chrom, start, stop, offsets).upper()
#	kmer_counts = template.copy()
#	for i in range(len(seq)-k+1):
#		kmer = seq[i:i+k]
#		if "N" in kmer: continue
#		kmer_counts[kmer] += 1
#		kmer_counts[rev_comp(kmer)] += 1
#	out = fields + [str(kmer_counts[x]) for x in kmers]
#	print "\t".join(out)
