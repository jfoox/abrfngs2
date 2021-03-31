#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00   # HH/MM/SS
#SBATCH --mem=32GB
#SBATCH --job-name=NGS2_DS_strelka2
#SBATCH --output=logs/NGS2_DS_strelka2.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

### get BAM
samplesheet=$1
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
reference_path=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)
knownsites_path=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz
rsync -av $knownsites_path* $TMPDIR
knownsites=$(basename $knownsites_path)
cd $TMPDIR

### run
STRELKA_INSTALL_PATH=/home/jfoox/programs/strelka-2.9.10.centos6_x86_64
${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam $bam \
    --referenceFasta $reference \
    --indelCandidates $knownsites \
    --runDir $PWD
runWorkflow.py -m local -j 4

### export
mv results ${sample}_strelka2
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs
rsync -avr ${sample}_strelka2 $outdir

echo [`date`] Finished job

