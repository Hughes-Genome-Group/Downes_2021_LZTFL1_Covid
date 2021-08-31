#!/bin/bash
#SBATCH --partition=batch
#SBATCH --job-name=snake_covid
#SBATCH --ntasks=1
#SBATCH --mem=12G
#SBATCH --time=02-00:00:00
#SBATCH --output=log_snake_covid_gwas.out
#SBATCH --error=log_snake_covid_gwas.err
#SBATCH --mail-user=ron.schwessinger@ndcls.ox.ac.uk
#SBATCH --mail-type=end,fail

cd /stopgap/covidgwas/rschwess/allelic_bias/wasp_pipe

module load python-cbrg/current
module load bowtie2
module load samtools

SLURM_ARGS="-p {cluster.partition} -n {cluster.ntasks} -t {cluster.time} -J {cluster.job-name} -o {cluster.output} -e {cluster.error}"

snakemake -j 8 -pr --cluster-config cluster.yaml --cluster "sbatch $SLURM_ARGS" --rerun-incomplete --latency-wait 120
