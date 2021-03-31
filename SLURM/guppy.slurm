#SBATCH --cpus-per-task=8
#SBATCH --time=168:00:00
#SBATCH --mem=64GB
#SBATCH --job-name=NGS2_guppy
#SBATCH --output=logs/NGS2_guppy.SAMPLE.%A_%a.log

source ~/.bashrc
echo [`date`] Started job

fast5dirfile=$1
fast5dir_path=$(cat $fast5dirfile | sed -n "${SLURM_ARRAY_TASK_ID}p")
rsync -avr $fast5dir_path $TMPDIR
fast5dir=$(basename $fast5dir_path)
outdir=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/ONT/outs
cd $TMPDIR

### run guppy
guppydir=/home/jfoox/programs/ont-guppy-cpu_4.0.11/bin
${guppydir}/guppy_basecaller \
  -i $fast5dir -s . \
  --flowcell FLO-MIN106 --kit SQK-LSK109 \
  --qscore_filtering -r --cpu_threads_per_caller 8 --pt_scaling -q 0

### export
mkdir -p ${outdir}/$fast5dir
rsync -av pass/* ${outdir}/$fast5dir

echo [`date`] Finished job


