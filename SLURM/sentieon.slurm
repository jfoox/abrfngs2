#!/bin/bash -l
#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=64GB
#SBATCH --job-name=NGS2_sen
#SBATCH --output=logs/NGS2_sen_%A_%a.log
source ~/.bashrc

echo [`date] Started job

### sentieon parameters and bwa
export SENTIEON_LICENSE=abclic01.med.cornell.edu:8990 
export SENTIEON_INSTALL_DIR=/home/jfoox/programs/sentieon-genomics-202010/bin/sentieon
export PATH=$PATH:/home/jfoox/programs/sentieon-genomics-202010/bin
export bwt_max_mem=48G

### metadata and refs
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/sentieon/outs
#rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
#rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/00-All.vcf.gz* $TMPDIR
reference=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
samplesheet=$1
fastqR1_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
fastqR1=$(basename $fastqR1_path)
sample=$(echo $fastqR1 | cut -d'.' -f1)

# account for single/paired end reactions
if [ -f ${fastqR1_path/.R1/.R2} ]; then fastqR2_path=${fastqR1_path/.R1/.R2}; else fastqR2_path=''; fi
if [ ! -z $fastqR2_path ]; then fastqR2=$(basename $fastqR2_path); fi
echo "/////////////////////////////////////////////////"
echo "This job's sample is: " $sample
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### initial alignment
if [ ! -f ${outdir}/${sample}.bam ] && [ ! -f ${outdir}/${sample}.md.bam ]; then
  rsync -avL $fastqR1_path $TMPDIR
  rsync -avL $fastqR2_path $TMPDIR
  echo "/////////////////////////////////////////////////"
  echo "*** Working on:"
  echo "*** $fastqR1"
  echo "*** $fastqR2"
  echo "/////////////////////////////////////////////////"
  
  sampleRG=$(echo $sample | cut -d'-' -f1)
  readgroup=@RG\\tID:${sample}\\tSM:${sampleRG}\\tPL:ILLUMINA
  sentieon bwa mem -M -R $readgroup -t 8 -K 10000000 $reference $fastqR1 $fastqR2 \
  |  sentieon util sort -r $reference -o ${sample}.bam -t 4 --sam2bam -i -
  rsync -av ${sample}.bam $outdir 
  samtools index -@3 ${sample}.bam
  rsync -av ${sample}.bam.bai $outdir 
  
  # calculate alignment metrics
  sentieon driver -t 4 -r $reference -i ${sample}.bam \
  --algo GCBias --summary GC_SUMMARY_TXT ${sample}_GC_METRIC.txt \
  --algo MeanQualityByCycle ${sample}_MQ_METRIC.txt \
  --algo QualDistribution ${sample}_QD_METRIC.txt \
  --algo InsertSizeMetricAlgo ${sample}_IS_METRIC.txt  \
  --algo AlignmentStat ${sample}_ALN_METRIC.txt
  rsync -av ${sample}_*.txt $outdir
fi

### mark duplicates
if [ ! -f ${outdir}/${sample}.md.bam ]; then
  if [ ! -f ${sample}.bam ]; then rsync -av ${outdir}/${sample}.bam* $TMPDIR; fi
  sentieon driver -t 4 -i ${sample}.bam --algo LocusCollector --fun score_info ${sample}_SCORE.gz
  sentieon driver -t 4 -i ${sample}.bam --algo Dedup --score_info ${sample}_SCORE.gz --metrics ${sample}_DEDUP_METRIC.txt ${sample}.md.bam
  rsync -av ${sample}.md.bam $outdir 
  samtools index -@3 ${sample}.md.bam
  rsync -av ${sample}.md.bam.bai $outdir 
  rsync -av ${sample}_DEDUP_METRIC.txt $outdir
fi

### variant calling
### (any sample with "-" has multiple input files (lanes, etc) and needs to be combined first)
if [[ ! $sample == *"-"* ]]; then
  if [ ! -f ${outdir}/${sample}.vcf ]; then
    if [ ! -f ${sample}.md.bam ]; then rsync -av ${outdir}/${sample}.md.bam* $TMPDIR; fi
    sentieon driver -t 4 -r $reference -i ${sample}.md.bam --algo QualCal ${sample}_RECAL_DATA.TABLE
    sentieon driver -t 4 -r $reference -i ${sample}.md.bam -q ${sample}_RECAL_DATA.TABLE --algo QualCal ${sample}_RECAL_DATA.TABLE.POST 
    sentieon driver -t 4 --algo QualCal --plot --before ${sample}_RECAL_DATA.TABLE --after ${sample}_RECAL_DATA.TABLE.POST ${sample}_RECAL_RESULT.CSV
    sentieon driver -t 4 -r $reference -i ${sample}.md.bam  -q ${sample}_RECAL_DATA.TABLE --algo Haplotyper --pcr_indel_model NONE -d 00-All.vcf.gz ${sample}.PCRindelNONE.vcf
    rsync -av ${sample}_*TABLE* $outdir
    rsync -av ${sample}.PCRindelNONE.vcf $outdir
  fi
else
  echo "$sample needs to be merged with its partners first using sambamba merge"
fi

echo [`date`] Finished job
