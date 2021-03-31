#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32GB
#SBATCH --job-name=NGS2_BactSen
#SBATCH --output=logs/NGS2_BactSen_%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### sentieon parameters and bwa
export SENTIEON_LICENSE=abclic01.med.cornell.edu:8990 
export SENTIEON_INSTALL_DIR=/home/jfoox/programs/sentieon-genomics-202010/bin/sentieon
export PATH=$PATH:/home/jfoox/programs/sentieon-genomics-202010/bin
export bwt_max_mem=12G

### metadata and refs
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/bacteria/sentieon/outs
samplesheet=$1
fastqR1_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d'@' -f1)
fastqR1=$(basename $fastqR1_path)
sample=$(echo $fastqR1 | cut -d'.' -f1)
if [[ $fastqR1 == *"R1"* ]]; then fastqR2_path=${fastqR1_path/R1/R2}; else fastqR2_path=''; fi
if [ ! -z $fastqR2_path ]; then fastqR2=$(basename $fastqR2_path); fi
reference_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d'@' -f2)
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)
single_or_pair=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p" | cut -d'@' -f3)

echo "/////////////////////////////////////////////////"
echo "This job's sample is:      " $sample
echo "This job's reference is:   " $reference
echo "This job's sample type is: " $single_or_pair
echo "/////////////////////////////////////////////////"
cd $TMPDIR

### initial alignment
if [ ! -f ${outdir}/${sample}.bam ] && [ ! -f ${outdir}/${sample}.md.bam ]; then
  rsync -avL $fastqR1_path $TMPDIR
  if [ ! -z $fastqR2_path ]; then rsync -avL $fastqR2_path $TMPDIR; fi
  echo "/////////////////////////////////////////////////"
  echo "*** Working on:"
  echo "*** $fastqR1"
  echo "*** $fastqR2"
  echo "/////////////////////////////////////////////////"
  
  if [[ $sample == "PGM"* ]] || [[ $sample == "S5"* ]]; then PLtag="IONTORRENT"; else PLtag="ILLUMINA"; fi
  readgroup=@RG\\tID:${sample}\\tSM:${sample}\\tPL:${PLtag}
  sentieon bwa mem -M -R $readgroup -t 4 -K 10000000 $reference $fastqR1 $fastqR2 \
  |  sentieon util sort -r $reference -o ${sample}.bam -t 4 --sam2bam -i -
  rsync -av ${sample}.bam $outdir 
  samtools index -@3 ${sample}.bam
  rsync -av ${sample}.bam.bai $outdir 
  
  # calculate metrics
  sentieon driver -t 4 -r $reference -i ${sample}.bam \
  --algo GCBias --summary GC_SUMMARY_TXT ${sample}_GC_METRIC.txt \
  --algo MeanQualityByCycle ${sample}_MQ_METRIC.txt \
  --algo QualDistribution ${sample}_QD_METRIC.txt \
  --algo InsertSizeMetricAlgo ${sample}_IS_METRIC.txt  \
  --algo AlignmentStat ${sample}_ALN_METRIC.txt
  rsync -av ${sample}_*.txt $outdir
fi

### mark duplicates (if paired end)
if [[ $single_or_pair == "paired" ]]; then
  if [ ! -f ${outdir}/${sample}.md.bam ]; then
    if [ ! -f ${sample}.bam ]; then rsync -av ${outdir}/${sample}.bam* $TMPDIR; fi
    sentieon driver -t 4 -i ${sample}.bam --algo LocusCollector --fun score_info ${sample}_SCORE.gz
    sentieon driver -t 4 -i ${sample}.bam --algo Dedup --score_info ${sample}_SCORE.gz --metrics ${sample}_DEDUP_METRIC.txt ${sample}.md.bam
    rsync -av ${sample}.md.bam $outdir 
    samtools index -@3 ${sample}.md.bam
    rsync -av ${sample}.md.bam.bai $outdir 
    rsync -av ${sample}_DEDUP_METRIC.txt $outdir
  fi
fi 

### variant calling
if [ ! -f ${outdir}/${sample}.vcf ]; then
  if [[ $single_or_pair == "paired" ]]; then bam=${sample}.md.bam; else bam=${sample}.bam; fi
  if [ ! -f $bam ]; then rsync -av ${outdir}/$bam* $TMPDIR; fi
  sentieon driver -t 4 -r $reference -i $bam --algo QualCal ${sample}_RECAL_DATA.TABLE
  sentieon driver -t 4 -r $reference -i $bam -q ${sample}_RECAL_DATA.TABLE --algo QualCal ${sample}_RECAL_DATA.TABLE.POST 
  sentieon driver -t 4 --algo QualCal --plot --before ${sample}_RECAL_DATA.TABLE --after ${sample}_RECAL_DATA.TABLE.POST ${sample}_RECAL_RESULT.CSV
  sentieon driver -t 4 -r $reference -i $bam  -q ${sample}_RECAL_DATA.TABLE --algo Haplotyper ${sample}.vcf
  rsync -av ${sample}_*TABLE* $outdir
  rsync -av ${sample}.vcf $outdir
fi

echo [`date`] "Finished job"

