#!/bin/bash -l
#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=16GB
#SBATCH --job-name=NGS2_mismatch-ucsc
#SBATCH --output=logs/NGS2_mismatch-ucsc.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### get BAM and BEDs
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avL $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/error/outs
ref=GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
cd $TMPDIR

### loop through GC bins
for num in $(seq 0.00 0.05 1.00); do
  if [[ $num == "0.00" ]]; then continue; fi
  GCmax=$num
  GCmin=$(bc <<< "$GCmax-0.04")
  GCmaxDecimal=$(echo $GCmax | cut -d'.' -f2)
  if [[ $GCmaxDecimal == "00" ]]; then GCmaxDecimal="100"; fi
  echo [`date`] " ///////////////////////////////// " $GCmax
  # need to create sub-BAM (or else stats will be on complete BAM)
  outbam=${sample}.0${GCmin}to${GCmax}.bam
  reformat.sh -in=$bam -out=$outbam mingc=$GCmin maxgc=$GCmax mappedonly=t primaryonly=t -ref=$ref -crashjunk=f -tossjunk=t -Xmx12g
  
  # now get statistics
  reformat.sh -in=$outbam mingc=$GCmin maxgc=$GCmax mappedonly=t primaryonly=t -ref=$ref -Xmx12g \
    ehist=ehist.${sample}.${GCmin}to${GCmaxDecimal}.txt \
    mhist=mhist.${sample}.${GCmin}to${GCmaxDecimal}.txt \
    indelhist=ihist.${sample}.${GCmin}to${GCmaxDecimal}.txt
  rsync -av ehist.*.txt ${outdir}/ehists
  rsync -av mhist.*.txt ${outdir}/mhists
  rsync -av ihist.*.txt ${outdir}/ihists
  
  # clear BAM to maintain scratch space
  rm $outbam
  
done

echo [`date`] 'Finished job'
