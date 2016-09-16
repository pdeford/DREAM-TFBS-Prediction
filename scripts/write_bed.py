#!/usr/bin/env python

"""
This script will read in a bed file that details the regions
being analyzed, and return a bed file containing overlapping
50bp regions in that range.

Usage: write_ped.py /cache/DREAM-tfbs-pred/annotations/regions/test_regions.blacklistfiltered.merged.bed 
"""


import sys

f = open(sys.argv[1])

chrom_sizes = {}
chrom_list = []

for line in f:
	fields = line.strip().split()
	chrom = fields[0]
	if chrom not in chrom_sizes: 
		chrom_sizes[chrom] = [int(fields[1]),None]
		chrom_list.append(chrom)
	chrom_sizes[chrom][1] = int(fields[2])

g = open(sys.argv[2], "wb")

count = 0
for chrom in chrom_list:
	bottom,top = chrom_sizes[chrom]
	for i in range(bottom-150,top+150,25):
		g.write("{0}\t{1}\t{2}\t{0}:{1}-{2}\n".format(chrom,i,i+50))

g.close()
