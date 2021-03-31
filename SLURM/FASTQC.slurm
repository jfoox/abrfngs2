#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=16GB
#SBATCH --job-name=NGS2_fastqc
#SBATCH --output=logs/NGS2_fastqc_%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### set up
samplesheet=$1
fastq_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $fastq_path $TMPDIR
fastq=$(basename $fastq_path)
sample=$(echo $fastq | cut -d'.' -f1)
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/fastqc/outs
cd $TMPDIR

### run
fastqc -t 4 $fastq
rsync -av *fastqc* $outdir

echo [`date`] Finished job
