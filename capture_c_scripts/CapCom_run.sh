#!/bin/bash
#$ -cwd
#$ -M ddownes
#$ -m eas
#$ -j n
#$ -N COVID_Comp

############################
###     Input Options    ###
############################

# Choice of comparison (either "2way" or "3way")              
analysis="3way"

# Directory containing 6/9 run files with structure described in README.md.
path="/path/to/CapSeqBasic/Results/"

# Desired output name of samples with the exact name of the respective input folders.
sample1="HUVEC"
sample1Directories="HUVEC_rep1,HUVEC_rep2,HUVEC_rep3"

sample2="Erythroid"
sample2Directories="Ery_rep1,Ery_rep2,Ery_rep3"

sample3="hESC" 
sample3Directories="ESC_rep1,ESC_rep2,ESC_rep3"

sample4="CD4a"
sample4Directories="CD4A_rep1,CD4A_rep2,CD4A_rep3"

sample5="CD4n"
sample5Directories="CD4N_rep1,CD4N_rep2,CD4N_rep3"

sample6="CD14" 
sample6Directories="CD14_rep2,CD14_rep3"

# Name for the analysis run: e.g. "GWAS_Promoter"
name="COVID"

# Genome (supports "hg18","hg19","mm9","mm10","dm3")
genome="hg38"

# CaptureC analyser version (Supports any version with output strcuture: F6_greenGraphs_combined_Samples_Version format) e.g. "CS5"
version="CS5"

# Path to file containing viewpoints and windowing parameters -
# Format: Viewpoint    Chr Frag_start  Frag_stop Exlusion_Start Exlusion_Stop Plot_Region_Start Plot_Region_Stop Bin_Size Window_size
parameters="/path/to/CapCom_parameters.txt"

# Path to where you would like the public UCSC hub to be made
public="/public/out/path/"

# Name of enzyme used: supports, "dpnII", "nlaIII", "hindIII"
enzyme="dpnII"

#====================================================#
### Core paths for files and scripts.
#====================================================#

#Path to cis3way shell.
CompareShell="/path/to/custom/cis3way.sh"

#Default annotation for plots is ref seq genes but any bed file can be inserted.
annotation="/path/to/annotationExamples/RefSeqGenes_${genome}.bed"

#Path to folder containing digested genomes - will search here for correct digest and create if not present.
digest="//restFrag_examples/"

#====================================================#
### Run Analysis
#====================================================#

# We are now here
rundir=$( pwd )
echo "Running ${CompareShell} in ${rundir}"

# Print run command
echo "${CompareShell} --analysis=${analysis} --path=${path} --sample1=${sample1} --sample2=${sample2} --sample3=${sample3} --sample4=${sample4} --sample5=${sample5} --sample6=${sample6} --directories=${sample1Directories},${sample2Directories},${sample3Directories},${sample4Directories},${sample5Directories},${sample6Directories} --name=${name} --genome=${genome} --version=${version} --parameters=${parameters} --annotation=${annotation} --public=${public} --restrictionenzyme=${enzyme} --frags=${digest}" 

# Run the command :
${CompareShell} --analysis=${analysis} --path=${path} --sample1=${sample1} --sample2=${sample2} --sample3=${sample3} --sample4=${sample4} --sample5=${sample5} --sample6=${sample6}  --directories=${sample1Directories},${sample2Directories},${sample3Directories},${sample4Directories},${sample5Directories},${sample6Directories} --name=${name} --genome=${genome} --version=${version} --parameters=${parameters} --annotation=${annotation} --public=${public} --restrictionenzyme=${enzyme} --frags=${digest}

echo "All done !"
echo


