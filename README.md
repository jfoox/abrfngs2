# ABRF NGS Phase II

## Introduction 
Analysis and figure generation code for the ABRF NGS Phase II Study on DNA-seq reproducibility. This repository includes scripts to run heavy lifting such as alignment and variant calling (SLURM), shell scripts to do post-processing calculations (bin), and R scripts used to create figures (Rmds).

## Reference materials
This study requires several resources in the `reference` directory. Included in this repo:
1. GRCh38 reference genome [`GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz`](ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz)
	* As recommended by Heng Li's [Which human reference genome to use?](https://lh3.github.io/2017/11/13/which-human-reference-genome-to-use)
2. RepeatMasker tracks from UCSC Table Browser
	* see [reference/ucsc/README_generate_UCSC_beds.txt](reference/ucsc/README_generate_UCSC_beds.txt)

Required:
1. Download dbSNP VCF [`GCF_000001405.25.gz`](https://ftp.ncbi.nlm.nih.gov/snp/latest_release/VCF/)
2. Download SDF reference for RTG analysis [`GRCh38.sdf.zip`](https://s3.amazonaws.com/rtg-datasets/references/GRCh38.sdf.zip)

## Primary analyses
SLURM scripts. You should run:

0. Quality Control with [`FASTQC.slurm`](SLURM/FASTQC.slurm)
1. Reference-based alignment
	* Cell line-specific alignment 
		* using Sentieon with [`sentieon.slurm`](SLURM/sentieon.slurm)
		* using GATK with [`GATK.slurm`](SLURM/GATK.slurm)
	* Bacteria-specific  alignment using Sentieon with [`sentieon_bacteria.slurm`](SLURM/sentieon_bacteria.slurm)
	* Get alignment statistics with [`BAMstats.slurm`](SLURM/BAMstats.slurm)
2. Downsample BAMs with [`downsampleBAMs.slurm`](SLURM/downsampleBAMs.slurm)
3. Variant calling
	* Sentieon-based Haplotyper with [`sentieonHaplotyper_on_downsampled.slurm`](SLURM/sentieonHaplotyper_on_downsampled.slurm)
	* Strelka2 with [`strelka2.slurm`](SLURM/strelka2.slurm)
	* DeepVariant with [`DeepVariant.slurm`](SLURM/DeepVariant.slurm)
	* GATK variant calling embedded in GATK script above
	