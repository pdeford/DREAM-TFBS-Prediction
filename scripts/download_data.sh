#!/usr/bin/env bash

cd /cache/DREAM-tfbs-pred/

# Download data from DREAM/Synapse
python ~/download_challenge_data.py /cache/DREAM-tfbs-pred/

# Create file structure
mkdir RNAseq

mkdir DNase
mkdir DNase/relaxed DNase/conservative DNase/signal

mkdir ChIP
mkdir ChIP/relaxed
mkdir ChIP/conservative ChIP/signal

mkdir annotations
mkdir annotations/labels annotations/regions

mkdir shape

# Move files to appropriate directories
mv ChIPseq.*signal* ChIP/signal/
mv ChIPseq.*conservative* ChIP/conservative/
mv ChIPseq.*relaxed* ChIP/relaxed/
mv gene_expression.* RNAseq/
mv *_regions*.bed annotations/regions/
mv *labels.tsv* annotations/labels/
mv DNASE.*conserv* DNase/conservative/
mv DNASE.*relax* DNase/relaxed/
mv DNASE.*signal* DNase/signal/
mv *_regions* annotations/regions/
gunzip hg19.genome.fa.gz 
gunzip annotations/labels/*

# Get and move additional files
wget "ftp://ftp.sanger.ac.uk/pub/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz"
gunzip gencode.v19.annotation.gtf.gz
mv gencode.v19.annotations* annotations

wget "ftp://rohslab.usc.edu/hg19/hg19.HelT.2nd.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.HelT.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.MGW.2nd.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.MGW.wig.bw"
#wget "ftp://rohslab.usc.edu/hg19/hg19.OC2.wig"
wget "ftp://rohslab.usc.edu/hg19/hg19.ProT.2nd.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.ProT.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.Roll.2nd.wig.bw"
wget "ftp://rohslab.usc.edu/hg19/hg19.Roll.wig.bw"

mv *.bw shape/

wget "https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes"
mv hg19.chrom.size annotations