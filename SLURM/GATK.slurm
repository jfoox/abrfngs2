#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=16GB
#SBATCH --job-name=NGS2_DS_GATK
#SBATCH --output=logs/NGS2_DS_GATK.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### GATK (based on GATK installation instructions: https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)
source activate gatk

### metadata and refs
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/00-All.vcf.gz* $TMPDIR
reference=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)

echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### variant calling
if [ ! -f ${outdir}/${sample}.gatk.vcf.gz ]; then
  if [ ! -f $bam ]; then rsync -av $bam_path* $TMPDIR; fi
  /home/jfoox/programs/gatk-4.1.9.0/gatk --java-options "-Xmx8g" \
    HaplotypeCaller --min-base-quality-score 10 --pcr-indel-model NONE --dbsnp 00-All.vcf.gz \
    -R $reference -I $bam -O ${sample}.gatk.vcf.gz
  rsync -av ${sample}.gatk.vcf.gz* $outdir
fi

echo [`date`] Finished job

