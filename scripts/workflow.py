#!/usr/bin/env python

"""
Brand new workflow for doing this stuff
"""

import cPickle as pickle
from multiprocessing import Pool
import random
import subprocess
import sys

import bx.bbi.bigwig_file
import numpy as np

###############################################################################

print >> sys.stderr, "Reading arguments"

TF = sys.argv[1]
cell = sys.argv[2]
data_dir = sys.argv[3] + "/"

training_labels_file = open(
        data_dir + "annotations/labels/{}.train.labels.tsv".format(TF)
    )

column_of_interest = [
    i for i,x in enumerate(training_labels_file.readline().split()) if x==cell
    ][0]

###############################################################################

def get_offsets():
    chrom_path =  data_dir + "hg19.genome.fa"

    offsets = {}

    with open(chrom_path) as f:
        while 1:
            line = f.readline()
            pos = f.tell()
            if line.startswith(">"): offsets[line.strip()[1:]] = pos
            elif line.strip()  ==  "": break

    return offsets

def lookup_sequence(chrom, start=None, end=None, offsets=None):
    chrom_path =  data_dir + "hg19.genome.fa"
    
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

nuc_index = dict(zip("ACGTN",range(5)))
def one_hot(sequence):
    oh = np.zeros([5, len(sequence)])
    for i,n in enumerate(sequence):
        oh[nuc_index[n], i] += 1
    return oh

def lookup_seq_structure(chrom, start, end):
    shape_dir = data_dir + "shape/"
    data = []
    for f_name in subprocess.check_output(["ls","{}".format(shape_dir)]).split():
        bwh = bx.bbi.bigwig_file.BigWigFile(open(shape_dir + "/" + f_name, "rb"))
        row = bwh.get_as_array(chrom, start, end)
        row[np.isnan(row)] = 0.0
        data.append(row)

    data = np.vstack(data)
    return np.ravel(data)

def norm_pdf(x,mu,var):
    l = (1./np.sqrt(2*np.pi*var))
    r = np.exp(-1*(x-mu)**2/(2*var))
    result = l*r
    #result = (1./np.sqrt(2*np.pi*var))*np.exp(-1*(x-mu)**2/(2*var))
    result += 10**-300
    return result

def bin(seq, binsize=50, func=max):
    out = []
    for i in range(0, len(seq), binsize):
        out.append(func(seq[i:i + binsize]))

###############################################################################
# Preload data

print >> sys.stderr, "Loading motifs"
# Motifs
pwm = pickle.load(open(data_dir + "out/{}_PWM.p".format(TF), 'rb'))
strum = pickle.load(open(data_dir + "out/{}_StruM.p".format(TF), 'rb'))

k = pwm.shape[1]
PWM = np.zeros([5, k])
PWM[:4,:] = pwm
PWM[4,:] += 0.25

p = len(subprocess.check_output(["ls", data_dir + "shape/"]).split())
pk = len(strum[0])/p


print >> sys.stderr, "Preparing DNase"
# DNase signal
dnase_signal_file = bx.bbi.bigwig_file.BigWigFile(
    open(data_dir + "DNase/signal/DNASE.%s.fc.signal.bigwig" % (cell), "rb")
    )

# DNase peaks
with open(
        data_dir + "DNase/conservative/DNASE.%s.conservative.narrowPeak" % (cell)
        ) as dnase_peak_file:
    peaks = {}
    for line in dnase_peak_file:
        fields = line.split()
        chrom, start, stop = fields[0], int(fields[1]), int(fields[2])
        if chrom not in peaks:
            peaks[chrom] = []
        peaks[chrom].append((start,stop))

###############################################################################

print >> sys.stderr, "Get offsets"
offsets = get_offsets()
frac_unbound = 0.05

Y = []#Y = []
data = []#data = []
#
for line in training_labels_file:#for line in training_labels_file:
#def wrapper(line):
    print line.strip()
    fields = line.split()
    bound = fields[column_of_interest]
    
    if bound == "A":
        continue#return #continue
    elif bound == "B":
        Y.append(1)#y = 1  #Y.append(1)
    elif random.random >= 1 - frac_unbound:
        Y.append(0)#y = 0  #Y.append(0)
    else:
        continue#return #continue

    chrom, start, stop = fields[0], int(fields[1]), int(fields[2])

    sequence = lookup_sequence(chrom, start - 150, stop + 150, offsets)
    one_hot_seq = one_hot(sequence)
    pwm_matches = []
    for i in range(len(sequence) - k + 1):
        kmer = one_hot_seq[:, i:(i + k)]
        kmer2 = kmer[:,::-1]
        kmer2[:4] = kmer2[3::-1]
        pwm_matches.append(
            max(np.sum(np.log(PWM*kmer)), np.sum(np.log(PWM*kmer2)))
            )

    struc_seq = lookup_seq_structure(chrom, start - 150, stop + 150)
    strum_matches = []
    for i in range(len(sequence) - pk + 1):
        kmer_struc = struc_seq[i*p:(i + pk)*p]
        kmer_struc2 = np.hstack(
            [kmer_struc[i:i+p] for i in range(0,len(kmer_struc),p)[::-1]]
            )
        try:
            score1 = np.sum(np.log2(norm_pdf(kmer_struc, strum[0], strum[1]**2)))
        except:
            print i
            quit()
        score2 = np.sum(np.log2(norm_pdf(kmer_struc2, strum[0], strum[1]**2)))
        strum_matches.append(
            max(score1, score2)
            )

    dnase_trace = dnase_signal_file.get_as_array(chrom, start - 150, stop + 150)
    dnase_trace[np.isnan(dnase_trace)] = 0.0
    
    pwm_feet = []
    for i, s in pwm_matches[50:-50 + k]:
        context = np.average([
                dnase_trace[i:(i + 50)], 
                dnase_trace[(i + 50 + k):(i + 100 + k)]
            ])
        center = np.average(dnase_trace[(i + 50):(i + 50 + k)])
        pwm_feet.append((context - center)*s)

    strum_feet = []
    for i, s in strum_matches[50:-50 + pk]:
        context = np.average([
                dnase_trace[i:(i + 50)], 
                dnase_trace[(i + 50 + pk):(i + 100 + pk)]
            ])
        center = np.average(dnase_trace[(i + 50):(i + 50 + pk)])
        pwm_feet.append((context - center)*s)


    for d_start, d_stop in peaks[chrom]:
        overlap = min(d_stop, stop) - max(d_start, start)
        if overlap > 0:
            dnase_overlap = overlap
            break
    else:
        dnase_overlap = 0

    #return 
    ([y, dnase_overlap, max(dnase_trace), max(pwm_matches), 
             max(strum_matches), max(pwm_feet), max(strum_feet)]
            + bin(dnase_trace, func=np.average) + bin(pwm_matches) 
            + bin(strum_matches) + bin(pwm_feet) + bin(strum_feet))

print >> sys.stderr, "Starting evaluation"
pool = Pool()
data = pool.map(wrapper, training_labels_file.readlines())
pool.close()
pool.join()

print >> sys.stderr, "Filtering"
Y = []
X = []
for i,row in enumerate(data):
    if row is not None:
        Y.append(row.pop(0))
        X.append(row)

Y = np.asarray(Y)
X = np.asarray(data)

print >> sys.stderr, "Done."

