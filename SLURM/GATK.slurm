#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=16GB

source ~/.bashrc
echo [`date`] Started job

### GATK (based on GATK installation instructions: https://gatk.broadinstitute.org/hc/en-us/articles/360035889851--How-to-Install-and-use-Conda-for-GATK4)
source activate gatk

### metadata and refs
rsync -av reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
rsync -av reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai $TMPDIR
rsync -av reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.dict $TMPDIR
rsync -av reference/00-All.vcf.gz* $TMPDIR
reference=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
outdir=$2

echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### variant calling
if [ ! -f ${outdir}/${sample}.gatk.vcf.gz ]; then
  if [ ! -f $bam ]; then rsync -av $bam_path* $TMPDIR; fi
  gatk --java-options "-Xmx8g" \
    HaplotypeCaller --min-base-quality-score 10 --pcr-indel-model NONE --dbsnp 00-All.vcf.gz \
    -R $reference -I $bam -O ${sample}.gatk.vcf.gz
  rsync -av ${sample}.gatk.vcf.gz* $outdir
fi

echo [`date`] Finished job

