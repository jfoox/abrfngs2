#SBATCH --time=168:00:00
#SBATCH --mem=32GB
#SBATCH --cpus-per-task=8

echo [`date`] Started job

### get BAM
samplesheet=$1
reference=reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
bam_path=$(cat $samplesheet | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -av $bam_path* $TMPDIR
bam=$(basename $bam_path)
sample=$(echo $bam | cut -d'.' -f1)
rsync -av deepvariant/deepvariant_1.1.0.sif $TMPDIR
cd $TMPDIR

### run
BIN_VERSION="1.1.0"
# CPU version:
singularity run -B /scratchLocal:/usr/lib/locale/ docker://google/deepvariant:"${BIN_VERSION}" \
  /opt/deepvariant/bin/run_deepvariant --model_type=WGS --ref=$reference --reads=$bam --output_vcf=${sample}.dv.vcf.gz \
  --intermediate_results_dir "intermediate_results_dir" --num_shards=8

### export
outdir=$2
rsync -av ${sample}* $outdir

echo [`date`] Finished job
