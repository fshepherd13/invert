# invert

### Pipeline for quantifying mRNA, cRNA and vRNA transcripts during influenza virus replication.

## Overview
Takes code from Thu Phan (https://github.umn.edu/phanx247/IAVKinetics) and works it into a reproducible pipeline using snakemake. Further details on the method can be found [here](https://pubmed.ncbi.nlm.nih.gov/33658346/).

## Usage
The pipeline requires the following dependencies:
* Snakemake and snakemake-wrapper-utils
* Trimmomatic
* STAR
* fastqc
* cufflinks
* samtools
* pandas

To install, use conda and the environment.yaml file within this repo:
`conda env create --file environment.yaml`

Activate the environment prior to running the pipeline:
`conda activate invert`

The user should only have to deal with adjusting the parameters within the `config.yaml` file. For detailed description, see the header "Config file" below.

To run the full pipeline, use:
`snakemake --cores 24` (<-- or however many cores you want Snakemake to utilize)

For further details on using Snakemake, see [here](https://snakemake.readthedocs.io/en/stable/)

## Config file
#### in_dir
This directory is the file path for where the raw fastq files are stored. 

#### samples
This is a user-specified list of sample names. In the case of UMGC at the U of M, each raw fastq file comes back with the format "sample name_R1_001.fastq.gz" for the forward read and "sample name_R2_001.fastq.gz" for the reverse read. Provide a comma-separated list that *only* contains the "sample name" part, i.e., do not include the "_R1_001.fastq.gz" or "_R2_001.fastq.gz". The pipeline will create output files that are named according to this sample list. 

#### Annotations and ref_genome
Although we are only interested in reads that are IAV specific, we do perform a STAR mapping step to map the trimmed reads to a concatenated genome of the host and IAV used in the experiment. In the subsequent step to quantify the mapped reads, only the IAV genome is used. These can be specified separately. See the header "preparing the reference files" for further details on preparing these files.

#### Genome_index
Directory that contains the index of the concatenated host (human, for human cell lines, or something else depending upon the cell line used in experiments) and IAV genome. This is produced by STAR. See "preparing the reference files" for details.

#### Parameters
These are user-specified parameters that tells certain programs how to be run. At the original writing of this, the only parameters here are included for the Trimmomatic step, which lists where the adapter sequences are located (the user has to change this file path! Adapter sequences are not included in this repo!) and the trailing option for Trimmomatic. 

## Preparing the reference files
There are several files that need to be prepared by the user:
* FASTA file of the IAV genome used in experiments
* FASTA file of the IAV genome concatenated to the genome of the host cell line used in experiments (i.e. IAV + human)
* GTF file of annotations for forward and reverse strands of IAV gene segments
* GTF file of the IAV annotations concatenated to the host cell line genome annotations (i.e. IAV + human)

First, generate fasta files for each of the 8 IAV segments. Save to a directory (I create a directory called `ref_files` to save my reference genomes and annotations). For example, your directory should contain the following files:
`Cal_HA.fa` 
`Cal_M.fa`
`Cal_NA.fa`
`Cal_NP.fa`
`Cal_NS.fa`
`Cal_PA.fa`
`Cal_PB1.fa`
`Cal_PB2.fa`
Concatenate them into a single file with a list of all 8 gene segments, using, for example:
`cat *.fa > IAV_Cal09.fa`
Clean up the directory by removing the individual gene segment files. See the example `IAV_Cal09.fa` file within `ref_files`.

Next, create an annotation file for the IAV strain that designates the forward and reverse segments. See the example `IAV_Cal09.gtf` file within `ref_files` and modify to fit your strain. Columns 1 and 2 should be edited to match the strain name, and columns 4-5 should be edited to match the positions within your segment (these positions must span the entire length of the gene, including the 5' and 3' UTRs). Transcript id's and gene id's should not be changed, these are set to attribute - sense strands to the viral genes (i.e. vPB2) while + sense are attributed to remaining rna (i.e. PB2).

Once you have a single fasta and gtf file for the IAV strain, you can use the included `compile_genome.srun` script. This script downloads the reference genome and gtf annotation for human (or whatever host cell line was used, this should be changed by the user!) from ensembl, concatenated with the IAV fasta and gtf files, and indexes them with STAR. Edit as needed to fit your situation. 