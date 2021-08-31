#!/bin/bash
set -e
set -u
set -o pipefail

# Filename HUVEC_rep1
# Align_pipe	Y
# MCC_splitter	Y
# MCCanal	Y
# MNase_Multi_limit	0
# MNase_multi	N
# Pipe2	Y
# bigwig_folder	/path/to/my_ref/
# blatN
# blat_command	
# bowtie	Y
# bowtie_genome_path	/databank/igenomes/Homo_sapiens/UCSC/hg38/Sequence/Bowtie2Index/genome
# clean_wig	N
# colour_file	/oligo_12_split_800bp.bed
# combine	0
# email	hidden@ndcls.ox.ac.uk
# flash	Y
# fq2fa	N
# genome	hg38
# gunzip	N
# gzip	N
# macs2	N
# master_folder	/path/to/00Scripts
# oligo_file	/path/to/oligo_12_split_800bp.fa
# postgzip	Y
# public_folder	/path/to/out
# public_folder_expt	/dpath/to/out
# public_url	https://sara.molbiol.ox.ac.uk
# qsub	1
# reference	COVID
# samtobam	Y
# sort	N
# track_hub	Y
# trim_galore	N
perl /path/to/MCC_pipe2.pl -p /path/to/HUVEC_rep3 -n HUVEC_rep3 -i MCC_pipe_config.txt -q
if wait; then echo "#########################################################################################
HUVEC_rep1 Pipe2 started
#########################################################################################">>log.txt;fi
