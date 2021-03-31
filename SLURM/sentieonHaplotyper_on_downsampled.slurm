#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=32GB
#SBATCH --job-name=NGS2_DS_sentieon
#SBATCH --output=logs/NGS2_DS_sentieon.SAMPLE.%A_%a.log

# run Sentieon Haplotyper on downsampled BAMs

source ~/.bashrc
echo [`date`] Started job

### sentieon parameters and bwa
export SENTIEON_LICENSE=abclic01.med.cornell.edu:8990 
export SENTIEON_INSTALL_DIR=/home/jfoox/programs/sentieon-genomics-202010/bin/sentieon
export PATH=$PATH:/home/jfoox/programs/sentieon-genomics-202010/bin
export bwt_max_mem=48G

rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/00-All.vcf.gz* $TMPDIR
rsync -av /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz* $TMPDIR
knownsites=resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs

### get BAM
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
reference_path=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)
cd $TMPDIR


### variant calling
sentieon driver -t 4 -r $reference -i $bam --algo QualCal -k $knownsites ${sample}_RECAL_DATA.TABLE
sentieon driver -t 4 -r $reference -i $bam -q ${sample}_RECAL_DATA.TABLE --algo QualCal -k $knownsites ${sample}_RECAL_DATA.TABLE.POST 
sentieon driver -t 4 --algo QualCal --plot --before ${sample}_RECAL_DATA.TABLE --after ${sample}_RECAL_DATA.TABLE.POST ${sample}_RECAL_RESULT.CSV
sentieon driver -t 4 -r $reference -i $bam  -q ${sample}_RECAL_DATA.TABLE --algo Haplotyper --pcr_indel_model NONE -d 00-All.vcf.gz ${sample}.sen.vcf
rsync -av ${sample}.sen.vcf $outdir

echo [`date`] Finished job
