#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=16GB

echo [`date`] Started job

### get BAM and BEDs
samplesheet=$1 # this is a text file containing the absolute path to one BAM per line
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
rsync -av reference/ucsc/ucsc.*.bed $TMPDIR
outdir=$2 # state which output directory you want
cd $TMPDIR

### run on whole BAM
samtools stats -@3 -t $bam > ${sample}.global.samstats

### loop through contexts
for bed in *.bed; do
  context=$(echo $bed | cut -d'.' -f2)
  echo [`date`] $context
  samtools stats -@3 -t $bed $bam > ${sample}.${context}.samstats
  rsync -av ${sample}.${context}.samstats $outdir
done

echo [`date`] Finished job