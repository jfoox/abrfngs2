# UCSC RepeatMasker regions retrieved from TableBrowser
# https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1047766977_eXZWt633c5adwac03fmCalNs7OfI
# (hg38, group: Repeats, track: Simple Repeats, table: simpleRepeat)
# saved as simpleRepeats_hg38.tsv
tail -n +2 simpleRepeats_hg38.tsv | grep -P '\t.$' | awk '{ print $2":"$3"-"$4 }' > simpleRepeats_homo.bed
tail -n +2 simpleRepeats_hg38.tsv | grep -vP '\t.$' | awk '{ print $2":"$3"-"$4 }' > simpleRepeats_STRs.bed
shuf -n 50000 simpleRepeats_homo.bed > simpleRepeats_homo.50k.bed
shuf -n 50000 simpleRepeats_STRs.bed > simpleRepeats_STRs.50k.bed

# get statistics within each region for each BAM 
function capture_errors {
  bam=$1
  sample=$(basename $bam | cut -d'.' -f1)
  regions=$2
  regiontype=$(echo $(basename $regions) | sed -e "s/simpleRepeat_//" | sed -e "s/.bed//")
  for i in `cat $regions`; do 
    if [[ $(samtools view $bam $i | head | wc -l) > 0 ]]; then
      samtools view -bS $bam $i \
      | samtools stats - \
      | grep 'error rate' \
      | cut -f3 
    else echo "NA"; fi 
  done > ${sample}.${regiontype}.out 2> /dev/null
}
export -f capture_errors
homo=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/simpleRepeats_homo.50k.bed
strs=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/simpleRepeats_STRs.50k.bed
ls /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs/*.bam | env_parallel -j 15 -v --lb "capture_errors {} $homo" &
ls /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs/*.bam | env_parallel -j 15 -v --lb "capture_errors {} $strs" &

# compile outputs into matrices
homo_tsv=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/simpleRepeats_homo.50k.tsv
strs_tsv=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/simpleRepeats_STRs.50k.tsv
outs_homo=$(ls *homo.50k.out | grep -v ALL)
outs_strs=$(ls *STRs.50k.out | grep -v ALL)
paste $homo_tsv $outs_homo \
  | sed -e "1s/^/bin,chrom,chromStart,chromEnd,name,period,copyNum,consensusSize,perMatch,perIndel,score,A,C,G,T,entropy,motif,$(echo $outs_homo | sed -e "s/.simpleRepeats//g" | sed -e "s/.homo.50k.out//g")\n/" | tr '\t' ',' | tr ' ' ',' \
  > ../../../tables/mismatch_homo_v2.csv
paste $strs_tsv $outs_strs  \
  | sed -e "1s/^/bin,chrom,chromStart,chromEnd,name,period,copyNum,consensusSize,perMatch,perIndel,score,A,C,G,T,entropy,motif,$(echo $outs_strs | sed -e "s/.simpleRepeats//g" | sed -e "s/.STRs.50k.out//g")\n/" | tr '\t' ',' | tr ' ' ','  \
   > ../../../tables/mismatch_strs_v2.csv
