# compile RTG outputs
# overall metrics (precision, sensitivity, F1) using best F value
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/rtg/outs
for i in $(ls | grep "rtg_" | grep -Ev "Alu|L1|L2|lowcom|LTR|satel|simple"); do 
  sample=$(echo $i | cut -d_ -f2-5)
  algo=$(echo $i | cut -d_ -f6)
  cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${algo}\t/g"
done | tr '\t' ',' | sed -e "1s/^/sample,caller,precision,sensitivity,fmeasure\n/" > ../../tables/rtg_summary_v2.csv

# clinvar
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/clinvar/outs
for i in *; do sample=$(echo $i | cut -d_ -f3-6); cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${algo}/g"; done | tr '\t' ',' | sed -e "1s/^/sample,precision,sensitivity,fmeasure\n/" > ../../tables/rtg_clinvar.csv
# OMIM
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/OMIM/outs
for i in *; do sample=$(echo $i | cut -d_ -f3-6); cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${algo}/g"; done | tr '\t' ',' | sed -e "1s/^/sample,precision,sensitivity,fmeasure\n/" > ../../tables/rtg_OMIM.csv
# exome
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/exome/outs
for i in *; do sample=$(echo $i | cut -d_ -f3-6); cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${algo}/g"; done | tr '\t' ',' | sed -e "1s/^/sample,precision,sensitivity,fmeasure\n/" > ../../tables/rtg_exome.csv


# ROC plots from SNP and INDEL tables
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/rtg/outs
for i in rtg_*; do 
  sample=$(echo $i | cut -d_ -f2-5)
  algo=$(echo $i | cut -d_ -f6)
  zgrep -v '#' ${i}/snp_roc.tsv.gz | sed -e "s/^/${sample}\t${algo}\t/g"
done | tr '\t' ',' | sed -e "1s/^/sample,caller,score,true_positives_baseline,false_positives,true_positives_call,false_negatives,precision,sensitivity,f_measure\n/" > ../../tables/rtg_ROC_SNP.csv

### UCSC RepeatMasker contexts
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/rtg/outs
for i in rtg_*_dv_*; do 
  sample=$(echo $i | cut -d_ -f2-5)
  context=$(echo $i | cut -d_ -f7)
  cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${context}\t/g"
done | tr '\t' ',' | sed -e "1s/^/sample,context,precision,sensitivity,fmeasure\n/" > ../../tables/rtg_ucsc.csv
# add long reads
cd /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/rtg/outs
for i in rtg_P*; do 
  sample=$(echo $i | cut -d_ -f2-5)
  context=$(echo $i | cut -d_ -f6)
  cat ${i}/summary.txt | sed -n "3p" | sed -r "s/ +/\t/g" | cut -f7-9 | sed -e "s/^/${sample}\t${context}\t/g"
done | tr '\t' ',' >> ../../tables/rtg_ucsc.csv


### get counts of TPs in each context for each genome
for bed in /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/ucsc/ucsc.*.bed; do 
  context=$(basename $bed | cut -d'.' -f2)
  for genome in Father Mother Son; do 
    if   [[ $genome == *"Mother"* ]]; then vcf=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/HG004_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
    elif [[ $genome == *"Father"* ]]; then vcf=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/HG003_GRCh38_1_22_v4.2.1_benchmark.vcf.gz
    elif [[ $genome == *"Son"*    ]]; then vcf=/athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/reference/HG002_GRCh38_1_22_v4.2.1_benchmark.vcf.gz; fi
    echo $context,$genome,SNPs,$(bedtools   intersect -a $vcf -b $bed | awk ' length($4) == length($5) ' | wc -l)
    echo $context,$genome,INDELs,$(bedtools intersect -a $vcf -b $bed | awk ' length($4) != length($5) ' | wc -l)
  done
done
