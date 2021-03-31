cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables

### mapping stats
echo "sample,readsTotal,readsProperlyPaired,readsDuplicated,readsMultimap,readsUnmap,errorRate,avgLength,avgQual" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_all_v2.csv
# short reads (manually edit Genapsys)
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/samstats/outs/full_BAMs
for i in $(ls *.samstats | grep -vE "PacBio|PromethION|MinION|Flongle|S5|Proton"); do
  sample=$(echo $i | cut -d'.' -f1)
  readsTotal=$(grep "raw total sequences" $i | cut -f3)
  readsPropPair=$(grep "reads properly paired:" $i | cut -f3)
  readsDup=$(grep "reads duplicated:" $i | cut -f3)
  readsMQ0=$(grep "reads MQ0:" $i | cut -f3)
  readsUnmap=$(grep "reads unmapped:" $i | cut -f3)
  errorRate=$(grep "error rate:" $i | cut -f3)
  avgLength=$(grep "average length:" $i | cut -f3)
  avgQual=$(grep "average quality:" $i| cut -f3)
  avgLength=$(grep "average length" $i | cut -f3)
  echo $sample,$readsTotal,$readsPropPair,$readsDup,$readsMQ0,$readsUnmap,$errorRate,$avgLength,$avgQual >> /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_all_v2.csv
done

# long reads
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/samstats/outs/full_BAMs
for i in $(ls *.samstats | grep -E "PacBio|PromethION|MinION|Flongle|S5|Proton"); do
  sample=$(echo $i | cut -d'.' -f1)
  readsTotal=$(grep "raw total sequences" $i | cut -f3)
  readsPropPair=$(grep "reads mapped:" $i | cut -f3)
  #readsDup=$(grep "reads duplicated:" $i | cut -f3)
  readsMQ0=$(grep "reads MQ0:" $i | cut -f3)
  readsUnmap=$(grep "reads unmapped:" $i | cut -f3)
  errorRate=$(grep "error rate:" $i | cut -f3)
  avgLength=$(grep "average length:" $i | cut -f3)
  avgQual=$(grep "average quality:" $i| cut -f3)
  echo $sample,$readsTotal,$readsPropPair,NA,$readsMQ0,$readsUnmap,$errorRate,$avgLength,$avgQual >> /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_all_v2.csv
done


# Bacterial

## short reads
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/samstats/outs/bacteria
echo "sample,readsTotal,readsProperlyPaired,readsDuplicated,readsMultimap,readsUnmap,errorRate,avgLength,avgQual" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_bacteria.csv
# paired end
for i in $(ls *.samstats | grep "MiSeq"); do
  sample=$(echo $i | cut -d'.' -f1)
  readsTotal=$(grep "raw total sequences" $i | cut -f3)
  readsPropPair=$(grep "reads properly paired:" $i | cut -f3)
  readsDup=$(grep "reads duplicated:" $i | cut -f3)
  readsMQ0=$(grep "reads MQ0:" $i | cut -f3)
  readsUnmap=$(grep "reads unmapped:" $i | cut -f3)
  errorRate=$(grep "error rate:" $i | cut -f3)
  avgLength=$(grep "average length:" $i | cut -f3)
  avgQual=$(grep "average quality:" $i| cut -f3)
  avgLength=$(grep "average length" $i | cut -f3)
  echo $sample,$readsTotal,$readsPropPair,$readsDup,$readsMQ0,$readsUnmap,$errorRate,$avgLength,$avgQual >> /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_bacteria.csv
done

## long reads
for i in $(ls *.samstats | grep -v "MiSeq"); do
  sample=$(echo $i | cut -d'.' -f1)
  readsTotal=$(grep "raw total sequences" $i | cut -f3)
  readsPropPair=$(grep "reads mapped:" $i | cut -f3)
  #readsDup=$(grep "reads duplicated:" $i | cut -f3)
  readsMQ0=$(grep "reads MQ0:" $i | cut -f3)
  readsUnmap=$(grep "reads unmapped:" $i | cut -f3)
  errorRate=$(grep "error rate:" $i | cut -f3)
  avgLength=$(grep "average length:" $i | cut -f3)
  avgQual=$(grep "average quality:" $i| cut -f3)
  echo $sample,$readsTotal,$readsPropPair,NA,$readsMQ0,$readsUnmap,$errorRate,$avgLength,$avgQual >> /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/mapping_stats_bacteria.csv
done


# ----------------------------------------------------------------- #


# Insert size distribution 
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/sentieon/outs
for i in *_IS_METRIC.txt; do 
  sampleIndiv=$(echo $i | sed -e "s/_IS_METRIC.txt//")
  sampleOverall=$(echo $sampleIndiv | cut -d'-' -f1)
  tail -n +6 $i | sed -e "s/^/${sampleOverall}\t${sampleIndiv}\t/g" | tr '\t' ',' 
done | sed -e "1s/^/sampleOverall,sampleIndiv,size,freq\n/" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/insert_size_all.csv

### GC Content (Normalized)
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/sentieon/outs
for i in *_GC_METRIC.txt; do 
  sampleIndiv=$(echo $i | sed -e "s/_GC_METRIC.txt//")
  sampleOverall=$(echo $sampleIndiv | cut -d'-' -f1)
  tail -n +3 $i | cut -f3,7 | head -101 | sed -e "s/^/${sampleOverall}\t${sampleIndiv}\t/g" | tr '\t' ','
done | sed -e "1s/^/sampleOverall,sampleIndiv,GC,normcov\n/" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/GC_normalizedCoverage_all.csv

#### GC Content (Normalized) (bacteria)
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/bacteria/sentieon/outs
for i in *_GC_METRIC.txt; do 
  sampleIndiv=$(echo $i | sed -e "s/_GC_METRIC.txt//")
  sampleOverall=$(echo $sampleIndiv | cut -d'-' -f1)
  tail -n +3 $i | cut -f3,5 | head -101 | sed -e "s/^/${sampleOverall}\t${sampleIndiv}\t/g" | tr '\t' ','
done | sed -e "1s/^/sampleOverall,sampleIndiv,GC,normcov\n/" > /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/tables/bacteria_GC_rawCoverage_all.csv


# Coverage metrics
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/coverage/outs
for context in Alu L1 L2 lowcomplexity LTR satellite simplerepeat; do 
  for i in *${context}.cov.bed; do 
    sample=$(echo $i | cut -d'.' -f1)
    head -10000 $i | sed -e "s/^/${context},${sample},/g"
  done
done | tr '\t' ',' | sed -e "1s/^/context,sample,chr,start,end,cov\n/" > ../../tables/ucsc-coverage_v7.csv

# Chromosomal Mapping
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/sentieon/outs
for i in $(ls *md.bam | grep -v '-'); do 
  sample=$(echo $i | cut -d'.' -f1)
  echo $sample | tr '\n' ','
  samtools idxstats $i | grep -v 'decoy' | cut -f3 | awk '{ sum += $0 } END { print sum }' | tr '\n' ','
  samtools idxstats $i | grep 'decoy' | cut -f3 | awk '{ sum += $0 } END { print sum }'
done > ../../tables/mapping_byChrType.csv
