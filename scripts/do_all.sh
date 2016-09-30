#!/usr/bin/env bash

#============================================================
# PREPARE THE ENVIRONMENT

# Read command line arguments
data_dir=$1 #'/cache/DREAM-tfbs-pred/'
out_dir=$2 #'/cache/DREAM-tfbs-pred/out/'

# Set the right environment
source activate DREAM-TFBS

#============================================================
# DOWNLOAD DATA

if [ ! -e $data_dir'/hg19.genome.fa' ]
	then
		echo "=================================================="
		echo "DOWNLOADING DATA"
		bash scripts/download_data.sh
fi

#============================================================
# PREPROCESS THE DATA

# Create windows to analyze
if [ ! -e 'average_regions.bed']
	then
		echo "=================================================="
		echo "EXTRACTING WINDOWS"
		python scripts/write_bed.py $data_dir'/annotations/regions/test_regions.blacklistfiltered.merged.bed' 'average_regions.bed'
fi

# Prepare the DNase files
if [ ! -e $out_dir'/HepG2_DNase.tsv']
	then
		echo "=================================================="
		echo "PROCESSING DNASE"
		python scripts/do_DNase.py $data_dir'/DNase/signal/' $out_dir
fi

bash prepare_DNase_overlap.sh $data_dir

# Learn all of the motifs
cut -f1 TF_cellTypes.txt | tail -n+2 > tfs.txt

echo "=================================================="
echo "EXTRACTING MOTIFS"
while read TF;
	do 
		if [ ! -d $out_dir'/'$TF'_meme' ]
			then
				cat $data_dir'/ChIP/conservative/ChIPseq.'*'.'$TF'.conservative.train.narrowPeak' | python scripts/learn_motifs.py $TF $out_dir $data_dir'/shape/'; 
		fi
	done < tfs.txt


# Score all motifs at all positions
echo "=================================================="
echo "SCORING MOTIFS"
while read TF;
	do
		if [ ! -e $out_dir'/'$TF'_PWM.wig' ]
			then
				if [ ! -e $out_dir'/'$TF'_PWM.bw' ]
					then
						# Score PWM at EVERY position genome wide
						python scripts/PWM_to_wig.py $TF $out_dir'/'$TF'_PWM.p' > $out_dir'/'$TF'_PWM.wig'

						# Convert PWM scores to bigWig
						wigToBigWig $out_dir'/'$TF'_PWM.wig' $data_dir'/annotations/hg19.chrom.sizes' $out_dir'/'$TF'_PWM.bw'
						rm $out_dir'/'$TF'_PWM.wig'

						# Score StruM at EVER position genome wide
						python scripts/StruM_to_wig.py $TF $out_dir'/'$TF'_StruM.p' $out_dir'/'$TF'_PWM.p' $data_dir'/shape/' > $out_dir'/'$TF'_StruM.wig'

						# Convert PWM scores to bigWig
						wigToBigWig $out_dir'/'$TF'_StruM.wig' $data_dir'/annotations/hg19.chrom.sizes'  $out_dir'/'$TF'_StruM.bw'
						rm $out_dir'/'$TF'_StruM.wig'

						# Extract data from these files:
						bigWigAverageOverBed -minMax $out_dir'/'$TF'_PWM.bw'   average_regions.bed $out_dir'/'$TF'_PWM.tsv'
						bigWigAverageOverBed -minMax $out_dir'/'$TF'_StruM.bw' average_regions.bed $out_dir'/'$TF'_StruM.tsv'
				fi
		fi
	done < tfs.txt


# Find the closest 5-genes to every 50bp window, and pull out their avg. FPKM in all cell types
if [ ! -e $out_dir'/RNA_vals.tsv' ]
	then
		echo "=================================================="
		echo "PROCESSING RNA"
		python scripts/RNA_closest.py $data_dir'/annotations/gencode.v19.annotation.gtf' $data_dir'/RNAseq/' $data_dir'/annotations/hg19.chrom.sizes' > $out_dir'/RNA_vals.tsv'
fi

# Count kmers in every position
if [ ! -e $out_dir'/kmers.tsv' ]
	then
		echo "=================================================="
		echo "COUNTING KMERS"
		python scripts/count_kmers.py average_regions.bed 5 >  $out_dir'/kmers.tsv'
fi


#============================================================
# CREATE THE ARRAYS FOR TRAINING THE MODELS
echo "=================================================="
echo "EXTRACTING TRAINING DATA"
python scripts/train_all.py $data_dir $out_dir

#============================================================
# GET PROBABILITIES FOR ALL MODELS
echo "=================================================="
echo "TRAINING MODELS AND GENERATING PREDICTIONS"
python scripts/eval_all.py $data_dir $out_dir
