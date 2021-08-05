# Multi-platform assessment of DNA sequencing performance in the ABRF Next-Generation Sequencing Study


## Introduction 
Analysis and figure generation code for the ABRF NGS Phase II Study on DNA-seq reproducibility. This repository includes scripts to run heavy lifting such as alignment and variant calling (SLURM), shell scripts to do post-processing calculations (bin), and R scripts used to create figures (Rmds).

## Requirements

* [BBMap](https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmap-guide/)
* [bedtools](https://bedtools.readthedocs.io/en/latest/index.html)
* [DeepVariant](https://github.com/google/deepvariant/blob/r1.1/docs/deepvariant-quick-start.md)
* [FASTQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/)
* [GATK4](https://gatk.broadinstitute.org/hc/en-us/articles/360036194592-Getting-started-with-GATK4)
* Guppy (see documentation on [Nanopore Community page](https://nanoporetech.com/))
* [minimap2](https://github.com/lh3/minimap2)
* [Picard](https://broadinstitute.github.io/picard/)
* [RTG Tools](https://github.com/RealTimeGenomics/rtg-tools)
* [samtools](https://github.com/samtools/samtools)
* [Sentieon](https://www.sentieon.com/)
* [Strelka2](https://github.com/Illumina/strelka)

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
Written in the form of SLURM scripts (see [SLURM Documentation](https://slurm.schedmd.com/documentation.html)), but the core code can easily be taken out and run directly.

0. **Quality Control** with [`FASTQC.slurm`](SLURM/FASTQC.slurm)
1. Reference-based **alignment**
	* Cell line-specific alignment 
		* using Sentieon with [`sentieon.slurm`](SLURM/sentieon.slurm)
		* using GATK with [`GATK.slurm`](SLURM/GATK.slurm)
	* Bacteria-specific  alignment using Sentieon with [`sentieon_bacteria.slurm`](SLURM/sentieon_bacteria.slurm)
	* Long read alignment with [`minimap2.slurm`](SLURM/minimap2.slurm)
	* Get alignment statistics with [`BAMstats.slurm`](SLURM/BAMstats.slurm)
	* See `tables-mapping.sh` and `tables-error.sh` below for creating tables with statistics based on alignments
2. **Downsample BAMs** with [`downsampleBAMs.slurm`](SLURM/downsampleBAMs.slurm)
3. **Variant calling**
	* Sentieon-based Haplotyper with [`sentieonHaplotyper_on_downsampled.slurm`](SLURM/sentieonHaplotyper_on_downsampled.slurm)
	* Strelka2 with [`strelka2.slurm`](SLURM/strelka2.slurm)
	* DeepVariant with [`DeepVariant.slurm`](SLURM/DeepVariant.slurm)
	* GATK variant calling embedded in GATK script above
4. **Other analyses**
	* Calculate mismatch rate in BAM with [`mismatchRate.slurm`](SLURM/mismatchRate.slurm)
	* Estimate variant detection sensitivity and precision with [`RTG.slurm`](SLURM/RTG.slurm)
	* Call FASTQs from ONT FAST5 data with [`guppy.slurm`](SLURM/guppy.slurm)
	
## Figure Generation
The code used to generate all figures (primary and Extended Data) are provided in the `Rmds` directory. (With the exception of Figure 1a, which was created in Adobe Illustrator.) `Rmds` includes:

1. Figure 1 (Depth of sequencing and mapping rate) with [`QCandMapping.R`](Rmds/QCandMapping.R)
	* Extended Data 1, 2 (see [`DecoyMapping.R`](Rmds/DecoyMapping.R))
2. Figure 2 (Genome Coverage) with [`Coverage.R`](Rmds/Coverage.R)
	* Extended Data 3, 10
3. Figure 3 (Mismatch rates) with [`Mismatch.R`](Rmds/Mismatch.R)
4. Figure 4 (Variant Detection) with [`Variants.R`](Rmds/Variants.R)
	* Extended Data 4, 5 (see [`VariantAlleles.R`](Rmds/VariantAlleles.R)), 6 (see [`MendelViolations.R`](Rmds/MendelViolations.R))
5. Figure 5 (Structural Variants) with [`SVs.R`](Rmds/SVs.R)
	* Extended Data 7, 8, 9
6. Figure 6 (Bacterial Sequencing) with [`Bacteria.R`](Rmds/Bacteria.R)

## Helper scripts
The `bin` directory contains python and shell scripts that enable primary analyses above. These include:

| Script | Function |
| ------ | -------- |
| **calculateMismatch.sh** | Calculate mismatch rates in homopolymer and STR contexts |
| **Mendel_upSetMatrixGen.py** | Create UpSet plots for Mendelian violations |
| **Mendel_violationsByType.py** | Parse outputs of VBT for Mendelian violations |
| **mendel.sh** | Run VBT Mendelian violations |
| **tables-error.sh** | Generate mismatch histograms via BBMap |
| **tables-mapping.sh** | Several functions to create tables with alignment statistics |
| **tables-variants.sh** | Several functions to create tables with variant detection statistics |
| **variantAllele_GTtoMatrix.py** | Convert genotype matrix to TSV for plotting |

## End notes

Please see XXX for publication.

The genome sequences in this study are available as EBV-immortalized B-lymphocyte cell lines (from Coriell) as well as from DNA (from Coriell and NIST). All data generated within this study from these genomes are publicly available on NCBI Sequence Read Archive (SRA) under the BioProject PRJNA646948, within accessions SRR12898279-12898354. 

You can cite our work as follows: [tk]
