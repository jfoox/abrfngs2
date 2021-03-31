#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=16GB
#SBATCH --job-name=NGS2_samstats-ucsc
#SBATCH --output=logs/NGS2_samstats-ucsc.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### get BAM and BEDs
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | sed -e "s/forbbmap_//" | cut -d'.' -f1)
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/ucsc/ucsc.*.bed $TMPDIR
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/samstats/outs
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
