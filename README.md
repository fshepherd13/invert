# invert

### Pipeline for quantifying mRNA, cRNA and vRNA transcripts during influenza virus replication.

## Overview
Takes code conceived of and originally written by Thu Phan (https://github.umn.edu/phanx247/IAVKinetics) and works it into a reproducible pipeline using snakemake. Further details on the method can be found [in her original publication](https://pubmed.ncbi.nlm.nih.gov/33658346/).

## Usage
The pipeline requires the following dependencies:
* Snakemake and snakemake-wrapper-utils
* Trimmomatic
* STAR
* fastqc
* cufflinks
* samtools
* pandas

This pipeline runs using snakemake. Each rule contains a conda environment with the software necessary to run it. To run the pipeline, you only need to install snakemake, i.e:
`conda env create -n snakemake -c bioconda snakemake`
`conda activate snakemake`

Each rule will create an isolated conda environment automatically as it runs.

To run the full pipeline, use:
`snakemake --cores 24` (<-- or however many cores you want Snakemake to utilize)

You can adjust several parameters within the `config.yaml` file to control how the pipeline runs. For detailed description, see the header "Config file" below.

For further details on using Snakemake, see [here](https://snakemake.readthedocs.io/en/stable/)

## Please read: using other IAV strains with pipeline
INVERT relies on genome positions within the IAV segments that may differ depending upon what strain is used. For mapping, you must make sure that your IAV reference annotation file matches your strain (see "Preparing the reference files" below.)

Within two scripts, there are samtools view commands that you will need to adjust if you are using something other than the PR8 strain contained in `ref_files/PR8.fa`. In `scripts/cmrna_count.sh`, lines 7-13 extract the 3' end of the IAV segments. Adjust these to match your strain name and 3' positions (in particular the end range, which should span the entirety of the 3' UTR). Within `scripts/cmrna_splicing.sh`, lines 10-22 look for the splice junctions of the M and NS segments. You must also adjust these to match the segment name and splice sites for your strain as these can differ. Apologies that this step is not yet generalized to work with multiple strains.

## Config file
#### in_dir
This directory is the file path for where the raw fastq files are stored. 

#### sample_file
This is the name of the csv file that contains the sample IDs (see below for file description).

#### Annotations and ref_genome
Although we are only interested in reads that are IAV specific, we do perform a STAR mapping step to map the trimmed reads to a concatenated genome of the host and IAV used in the experiment. See the header "preparing the reference files" for further details on preparing these files.

#### Genome_index
Directory that contains the index of the concatenated host (human, for human cell lines, or something else depending upon the cell line used in experiments) and IAV genome. This is produced by STAR. See "preparing the reference files" for details.

#### Parameters
These are user-specified parameters that tells certain programs how to be run. At the original writing of this, the only parameters here are included for the Trimmomatic step, which lists where the adapter sequences are located (adapter sequences are included in this repo under `ref_files` but please ensure these are correct for your sequencing run!) and the trailing option for Trimmomatic. 

## `samples.csv` file
A user-specified list of sample names. In the case of UMGC at the U of M, each raw fastq file comes back with the format "sample-name_R1_001.fastq.gz" for the forward read and "sample name_R2_001.fastq.gz" for the reverse read. Provide a comma-separated file that *only* contains the "sample-name" part, i.e., do not include the "_R1_001.fastq.gz" or "_R2_001.fastq.gz". The pipeline will create output files that are named according to this sample list. 

## Preparing the reference files
There are several files that need to be prepared by the user:
* FASTA file of the IAV genome used in experiments
* FASTA file of the IAV genome concatenated to the genome of the host cell line used in experiments (i.e. IAV + MDCK)
* GTF file of annotations for forward and reverse strands of IAV gene segments
* GTF file of the IAV annotations concatenated to the host cell line genome annotations
* Index (generated by STAR) of the concatenated host and IAV genome

First, generate fasta files for each of the 8 IAV segments. The PR8 sequence used in Phan et al., 2021 is located in the `ref_files` directory. If you are re-running invert for a different IAV strain, you will need to create your own fasta file of IAV segments. For example, your directory may contain the following files (shown for a different strain, for example Cal09):
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
Clean up the directory by removing the individual gene segment files.

Next, create an annotation file for the IAV strain that designates the forward and reverse segments. An example annotation file for PR8 used by Thu Phan is found within `ref_files/`. You must modify this file to fit your strain and its ORFs! Columns 1 and 2 should be edited to match the strain name, and columns 4-5 should be edited to match the positions within your segment (these positions must span the entire length of the gene, including the 5' and 3' UTRs). Transcript id's and gene id's should not be changed, these are set to attribute negative sense strands to the viral genes (i.e. vPB2) while positive sense are attributed to mRNA and cRNA (i.e. PB2).

Once you have a single fasta and gtf file for the IAV strain used in your experiments, you can refer to the included `workflow/scripts/compile_genome.srun` script. This script guides the user through downloading the reference genome and gtf annotation for the host cell line (shown in script for Canis familiaris, for MDCK cells) from ensembl, concatenated with the IAV fasta and gtf files, and indexes them with STAR. Edit as needed to fit your situation.