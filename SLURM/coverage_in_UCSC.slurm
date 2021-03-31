#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=16GB
#SBATCH --job-name=NGS2_coverage-ucsc
#SBATCH --output=logs/NGS2_coverage-ucsc.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### get BAM and BEDs
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
#rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.genome $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/ucsc/ucsc.*.bed $TMPDIR
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/coverage/outs
cd $TMPDIR

### subset BAM to primary alignments only
echo [`date`] Creating subset of BAM with only primary, high quality mappings ...
samtools view -b -F 3844 -q 10 $bam > ${sample}.primary.bam

### loop through contexts
for bed in *.bed; do
  context=$(echo $bed | cut -d'.' -f2)
  echo [`date`] $context
  bedtools coverage -a $bed -b ${sample}.primary.bam -mean -sorted > ${sample}.${context}.cov.bed
  rsync -av ${sample}.${context}.cov.bed $outdir
done

echo [`date`] Finished job
