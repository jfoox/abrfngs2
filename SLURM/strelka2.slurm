#!/bin/bash -l
#SBATCH --cpus-per-task=4
#SBATCH --time=168:00:00
#SBATCH --mem=32GB

echo [`date`] Started job

### get BAM
samplesheet=$1 # this is a text file containing the absolute path to one BAM per line
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
reference_path=reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
rsync -av $reference_path* $TMPDIR
reference=$(basename $reference_path)
knownsites_path=reference/resources_broad_hg38_v0_Homo_sapiens_assembly38.known_indels.filtered.vcf.gz
rsync -av $knownsites_path* $TMPDIR
knownsites=$(basename $knownsites_path)
cd $TMPDIR

### run
# see https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/quickStart.md for installation
STRELKA_INSTALL_PATH=/path/to/strelka-2.9.10.centos6_x86_64
${STRELKA_INSTALL_PATH}/bin/configureStrelkaGermlineWorkflow.py \
    --bam $bam \
    --referenceFasta $reference \
    --indelCandidates $knownsites \
    --runDir $PWD
runWorkflow.py -m local -j 4

### export
mv results ${sample}_strelka2
outdir=$2 # state which output directory you want
rsync -avr ${sample}_strelka2 $outdir

echo [`date`] Finished job