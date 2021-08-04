#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32GB

# this script can be modified to be run in CLINVAR genes, OMIM genes, exomic regions, etc.
# below is code for analysis in UCSC RepeatMasker regions

echo [`date`] Started job

samplefile=$1
outdir=$2
sample=$(cat $samplefile | sed -n "${SLURM_ARRAY_TASK_ID}p")
dsdir=downsample/outs

# setup
rsync -avr reference/GCA_000001405.15_SDF $TMPDIR
reference=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna

# get proper reference
if [[ $sample == *"Mother"* ]]; then 
  truthfile_path=reference/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  bedfile_path=reference/HG004_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
elif [[ $sample == *"Father"* ]]; then
  truthfile_path=reference/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  bedfile_path=reference/HG003_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
elif [[ $sample == *"Son"* ]]; then
  truthfile_path=reference/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
  bedfile_path=reference/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed
fi
rsync -av $truthfile_path* $TMPDIR
truthfile=$(basename $truthfile_path)
rsync -av $bedfile_path $TMPDIR
bedfile=$(basename $bedfile_path)
cd $TMPDIR

readtype=$2
if [[ $readtype == "short" ]]; then
  queryfile=${sample}.dv.vcf.gz; else
  queryfile=${sample}.sorted.ds25.noMono.vcf.gz
fi
rsync -avL ${dsdir}/${queryfile}* $TMPDIR

# run through all UCSC RepeatMasker contexts
rsync -av reference/ucsc/*.bed $TMPDIR
for ucscbed in ucsc.*.bed; do
  context=$(echo $ucscbed | cut -d'.' -f2)
  echo [`date`] $context
  bedtools intersect -a $ucscbed -b $bedfile > ucscbed_subset.bed
  
  # run RTG
  rtg RTG_MEM=24G vcfeval -b $truthfile -c $queryfile -e ucscbed_subset.bed \
    --bed-regions ucscbed_subset.bed -t GCA_000001405.15_SDF -T 4 -o rtg_${sample}_${context}
  # export
  rsync -avr rtg_${sample}_${context} $outdir
done

echo [`date`] Finished job
