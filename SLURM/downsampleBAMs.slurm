#SBATCH --cpus-per-task=1
#SBATCH --time=168:00:00
#SBATCH --mem=16GB

echo [`date`] Started job

### metadata and refs
bamdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/sentieon/outs
rsync -av reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna $TMPDIR
rsync -av reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.fai $TMPDIR
# autosomal means calculated using mosdepth with -n flag and canonical chromosome coverage averaged
rsync -avL all_autosomal_means.txt $TMPDIR
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

cov=$(grep $sample all_autosomal_means.txt | cut -d',' -f2)
avg=25
if [[ $cov < 25 ]]; then
  echo "Can't downsample because less than 25x: " $sample " /// " $cov
  exit
fi

## downsample
rsync -av $bam_path* $TMPDIR
prob=$(echo "scale=2; $avg / $cov" | bc)
picard=/path/to/picard.jar
java -jar $picard DownsampleSam \
    I=$bam \
    O=${bam/.bam/}.ds${avg}.bam \
    PROBABILITY=$prob \
    CREATE_INDEX=true \
    VALIDATION_STRINGENCY=LENIENT

# send back
mv ${bam/.bam/}.ds${avg}.bai ${bam/.bam/}.ds${avg}.bam.bai
rsync -av ${bam/.bam/}.ds${avg}.bam* $outdir

echo [`date`] Finished job
