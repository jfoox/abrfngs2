#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=NGS2_BactSen_ONT
#SBATCH --output=logs/NGS2_BactSen_ONT_%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### metadata and refs
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/bacteria/minimap2/outs
samplesheet=$1
fastqR1_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d'@' -f1)
fastqR1=$(basename $fastqR1_path)
sample=$(echo $fastqR1 | cut -d'.' -f1)
reference_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d'@' -f2)
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)

echo "/////////////////////////////////////////////////"
echo "This job's sample is:      " $sample
echo "This job's reference is:   " $reference
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### initial alignment
if [ ! -f ${outdir}/${sample}.bam ]; then
  rsync -avL $fastqR1_path $TMPDIR
  echo "/////////////////////////////////////////////////"
  echo "*** Working on:"
  echo "*** $fastqR1"
  echo "/////////////////////////////////////////////////"
  minimap2 --MD -a -t 8 -x map-ont -R @RG\tID:${sample}\tSM:${sample}\tPL:ONT $reference $fastq | sambamba view -S -f bam -t 8 -o ${sample}.bam /dev/stdin
  rsync -av ${sample}.bam $outdir 
  samtools index -@3 ${sample}.bam
  rsync -av ${sample}.bam.bai $outdir 
fi
echo [`date`] Finished job


