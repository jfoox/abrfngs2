#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32GB

# run Sentieon Haplotyper on downsampled BAMs

echo [`date`] Started job

### sentieon parameters and bwa
export SENTIEON_LICENSE=/path/to/license/server
export SENTIEON_INSTALL_DIR=/path/to/sentieon/bin/sentieon
export PATH=$PATH:/path/to/sentieon/bin
export bwt_max_mem=48G

rsync -av reference/00-All.vcf.gz* $TMPDIR
rsync -av reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz* $TMPDIR
knownsites=resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz

### get BAM
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
reference_path=reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)
outdir=$2 # state which output directory you want
cd $TMPDIR


### variant calling
sentieon driver -t 4 -r $reference -i $bam --algo QualCal -k $knownsites ${sample}_RECAL_DATA.TABLE
sentieon driver -t 4 -r $reference -i $bam -q ${sample}_RECAL_DATA.TABLE --algo QualCal -k $knownsites ${sample}_RECAL_DATA.TABLE.POST 
sentieon driver -t 4 --algo QualCal --plot --before ${sample}_RECAL_DATA.TABLE --after ${sample}_RECAL_DATA.TABLE.POST ${sample}_RECAL_RESULT.CSV
sentieon driver -t 4 -r $reference -i $bam  -q ${sample}_RECAL_DATA.TABLE --algo Haplotyper --pcr_indel_model NONE -d 00-All.vcf.gz ${sample}.sen.vcf
rsync -av ${sample}.sen.vcf $outdir

echo [`date`] Finished job