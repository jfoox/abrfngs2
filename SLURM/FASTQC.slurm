#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=16GB

echo [`date`] Started job

### set up
samplesheet=$1
fastq_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $fastq_path $TMPDIR
fastq=$(basename $fastq_path)
sample=$(echo $fastq | cut -d'.' -f1)
outdir=$2
cd $TMPDIR

### run
fastqc -t 4 $fastq
rsync -av *fastqc* $outdir

echo [`date`] Finished job
