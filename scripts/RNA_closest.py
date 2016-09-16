#!/usr/bin/env python

"""
Usage: script.py annotation.gtf data/RNAseq/ hg19.chrom.size > output.tsv
scripts/RNA_closest.py /cache/DREAM-tfbs-pred/annotations/gencode.v19.annotation.gtf /cache/DREAM-tfbs-pred/RNAseq/ /cache/DREAM-tfbs-pred/annotations/hg19.chrom.sizes > ooooout/RNA_vals.tsv 
"""

import sys
import numpy as np
import subprocess
from multiprocessing import Pool

############################################################

gtf_file = sys.argv[1] #/cache/DREAM-tfbs-pred/annotations/gencode.v19.annotation.gtf
RNA_path = sys.argv[2] #/cache/DREAM-tfbs-pred/RNAseq/
sizeFile = open(sys.argv[3]) #/cache/DREAM-tfbs-pred/annotations/hg19.chrom.sizes

# LOAD RNA SEQ DATA
print >> sys.stderr, "Loading RNA seq"
print >> sys.stderr, "...genes"
genes = {}
with open(gtf_file) as GTF:
	for line in GTF:
		if not line.startswith("#") and "transcript" in line:
			fields = line.split()
			genes[fields[9][1:-2]] = [fields[0], (int(fields[3])+int(fields[4]))/2] # genes[gene] = [chrom, midpoint]

RNA_files  = subprocess.check_output(["ls", RNA_path]).split()
cell_types = set([name.split(".")[1] for name in RNA_files])

print >> sys.stderr, "...expression"
for cell in cell_types:
	tmp = {}
	for name in RNA_files:
		if cell in name:
			with open(RNA_path+name) as RNA:
				RNA.readline()
				for line in RNA:
					fields = line.split()
					g = fields[0]
					FPKM = float(fields[-1])
					if g not in tmp:
						tmp[g] = []
					tmp[g].append(FPKM)
	for gene in tmp: genes[gene].append(np.average(tmp[gene]))
del tmp

# # Load expression values
# print "...test expression"
# with open(RNA_path + "/gene_expression.HepG2.biorep1.tsv") as RNA:
# 	RNA.readline()
# 	for line in RNA:
# 		fields = line.split()
# 		genes[fields[0]].append(float(fields[-1]))

# Extract positions and expression
print >> sys.stderr, "Reformatting"
RNA_positions,RNA_data = {},{}
gene_order = {}
for gene in genes.keys():
	if len(genes[gene]) < 3: 
		del genes[gene]
	else:
		data = genes[gene]
		chrom,mid,FPKM = data[0],data[1],data[2:]
		RNA_data[gene] = FPKM
		if chrom not in RNA_positions:
			RNA_positions[chrom],gene_order[chrom] = [],[]
		RNA_positions[chrom].append(mid)
		gene_order[chrom].append(gene)
for chrom in RNA_data: RNA_data[chrom] = np.asarray(RNA_data[chrom])
for chrom in RNA_positions: RNA_positions[chrom] = np.asarray(RNA_positions[chrom])
for chrom in gene_order: gene_order[chrom] = np.asarray(gene_order[chrom])
del genes

############################################################

def closestGenes(chrom,start,stop):
	sameChrom = RNA_positions[chrom]
	mid = (start+stop)/2
	distances = np.absolute(sameChrom-mid)
	winners = np.argsort(distances)[:5]
	return gene_order[chrom][winners]

############################################################

print >> sys.stderr, "Loading sizes"
sizes = {}
for line in sizeFile:
	fields = line.strip().split()
	sizes[fields[0]] = int(fields[1])

############################################################
print >> sys.stderr, "Processing"
chroms = ["chr{}".format(x) for x in range(1,23)+["X"]]

#print "#chrom\tstart\tstop\tgene1\tgene2\tgene3\tgene4\tgene5"

def wrapper(i):
	genes = list(closestGenes(chrom,i,i+50))
	values = [str(x) for x in np.hstack([RNA_data[gene] for gene in genes])]
	newLine = [chrom,str(i),str(i+50)] + values
	newLine = "\t".join(newLine)
	return newLine

for chrom in chroms:
	print >> sys.stderr, "...%s" % chrom
	#for i in range(0,sizes[chrom],50):
	#	genes = list(closestGenes(chrom,i,i+50))
	#	values = [str(x) for x in np.hstack([RNA_data[gene] for gene in genes])]
	#	newLine = [chrom,str(i),str(i+50)] + values
	#	newLine = "\t".join(newLine)
	#	#print newLine
	pool = Pool()
	newlines = pool.map(wrapper,range(0,sizes[chrom],50))
	print "\n".join(newlines)
