cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/mendelian_violation
git clone https://github.com/sbg/VBT-TrioAnalysis.git
cd VBT-TrioAnalysis
make all

function mendelian_analysis {
  vcfpath=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs
  ref=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
  mother=$1
  father=${mother/Mother/Father}
  son=${mother/Mother/Son}
  output=$(echo $(basename $mother) | cut -d'.' -f1)
  ./vbt mendelian \
    -ref $ref \
    -mother $mother \
    -father $father \
    -child $son \
    -outDir /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/mendelian_violation/outs \
    -out-prefix $output \
    --output-violation-regions \
    --autosome-only \
    -bed /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/HG002_GRCh38_1_22_v4.2.1_benchmark_noinconsistent.bed \
    -thread-count 2
}
export -f mendelian_analysis
  
mothers=$(ls /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs/*Mother_REP01.dv.vcf.gz)
echo $mothers | tr ' ' '\n' | env_parallel --lb -v -j 4 "mendelian_analysis {}"

# concatenate output
for i in *_tab_delim_detailed_log.tsv; do
  sample=$(echo $i | sed -e "s/_tab_delim_detailed_log.tsv//")
  rate=$(grep -m1 'Violation Rate' ${i/_tab_delim_detailed_log.tsv/_DetailedLogs.txt} | cut -d':' -f2)
  cat <(cat $i | cut -f3 | sed -e "s/VIOLATION/${sample}/") <(echo $rate) > col_${sample}
done
paste <(cut -f 1 HiSeq2500_LAB01_tab_delim_detailed_log.tsv | sed -e '$aViolationRate') col_* | less -S > ../../tables_mendelian_table.tsv

# remove unnecessary part of filenames
for i in *_Mother_REP01_*; do mv -i $i ${i/_Mother_REP01/}; done

# custom script for creating BED files based on region type
for i in *_trio.vcf; do python get_violations_by_type.py $i; done

# create UpSet plot matrices
python run_upSet_matrix_generator.py
