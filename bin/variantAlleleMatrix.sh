# merge variant call sets for short (DeepVariant) and long (clair2) reads
bcftools merge --threads 8 $(ls /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs/*dv.vcf.gz) $(ls /athena/masonlab/scratch/projects/abrf_ngs_phaseii/revision/downsample/outs/*.clair.noMono.vcf.gz) > all_WGS.vcf

# reduce to chr1
cat <(grep '#' all_WGS.vcf) <(grep -P 'chr1\t' all_WGS.vcf) > all_WGS.chr1.vcf

# convert to matrix
bcftools query -i 'MIN(DP)>9' -f '%CHROM\t%POS\t%REF\t%ALT\t[%GT\t]\n' all_WGS.chr1.vcf > all_WGS.chr1.GT
header=$(grep -m1 '#CHROM' all_WGS.chr1.vcf | cut -f10- | tr '\t' '@' | sed -e "s/^/Chr@Pos@Ref@Alt@/")
cat all_WGS.chr1.GT | sed -e "1s/^/$header\n/" | tr '@' '\t' > all_WGS.chr1.GT2
mv all_WGS.chr1.GT2 all_WGS.chr1.GT
