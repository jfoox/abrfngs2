# ABRF NGS Phase II
Analysis and figure generation code for the ABRF NGS Phase II Study on DNA-seq reproducibility. This repository includes scripts to run heavy lifting such as alignment and variant calling (SLURM), shell scripts to do post-processing calculations (bin), and R scripts used to create figures (Rmds).

### reference materials

This study uses common reference materials including:
* GRCh38 reference genome:  ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz 
** As recommended by Heng Li https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use
* UCSC RepeatMasker tracks from Table Browser
** see reference/ucsc/README_generate_UCSC_beds.txt
 