#!/bin/bash

#$ -cwd
#$ -M rschwess
#$ -m eas
#$ -j y
#$ -N cov_dephaem_4k

# Run as:
# qsub <pathtoscript.sh>

### START ###
 
date 
 
OUTDIR=/t1-data/user/hugheslab/rschwess/other_projects/covid_damien

cd ${OUTDIR}

#module load python-deeplearning/201702
module load bedtools

python /t1-data/user/hugheslab/rschwess/scripts/machine_learning/epigenome_nets/deepHaem/run_deploy_net.py --dlmodel deepHaemWindow \
	--batch_size 3 \
	--out_dir ./deepHaem_4k_haplo_out \
	--name_tag deepHaem_4k_haplo_predict --do damage_and_scores \
	--input rs17713054_haplotype_for_deephaem.vcf \
	--model /t1-data/user/hugheslab/rschwess/machine_learning/deepHaem/training_runs_archive/basenji_compendium_runs/chrom_full/train_dhw_plain1000_chrom_full_full_run_28022020/best_checkpoint-794328 \
	--genome /t1-data/user/hugheslab/rschwess/database/whole_genome_fasta/hg38.fa \
	--bp_context 1000 \
	--num_classes 4384 \
	--run_on cpu

echo "Finished ..."

date
