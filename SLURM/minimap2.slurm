#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=NGS2_PacBio
#SBATCH --output=logs/NGS2_PacBio_%A_%a.log
source ~/.bashrc

echo [`date`] Started job

### metadata and refs
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/minimap2/outs
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai $TMPDIR
ref=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
samplesheet=$1
fastqR1_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $fastqR1_path $TMPDIR
fastqR1=$(basename $fastqR1_path)
library=$(echo $fastqR1_path | rev | cut -d'/' -f2 | rev)
sample=${library}_${SLURM_ARRAY_TASK_ID}
echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"
cd $TMPDIR

minimap2 -a -t 8 -x map-pb --MD -H -R @RG\tID:${library}_${SLURM_ARRAY_TASK_ID}\tSM:${library}\tPL:PacBio $ref $fastqR1 | sambamba view -S -f bam -t 8 -o ${sample}.bam /dev/stdin
sambamba sort -t 8 ${sample}.bam

rsync -av ${sample}.sorted.bam* $outdir

echo [`date`] 'Finished job'
