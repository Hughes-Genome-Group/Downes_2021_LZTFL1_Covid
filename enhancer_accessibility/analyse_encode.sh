#!/bin/bash
#$ -cwd
#$ -M ddownes
#$ -m eas
#$ -j n
#$ -N ENCODE_corr




bed="promoter_correl.bed"

module load deeptools
module list

names=""
list=""

for file in *bw
do
       
        shorter=$( echo $file | sed 's/.bw//' )
        names="${names} ${shorter}"
        list="${list} ${file}"
done

echo "${list}"
echo "${names}"


multiBigwigSummary BED-file -b ${list} -o ENCODE_95_Coverage_correl.npz --BED ${bed} --labels ${names} --outRawCounts ENCODE_95_Coverage_correl.tab

exit

