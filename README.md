# README

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4741499.svg)](https://doi.org/10.5281/zenodo.4741499)  

**This pipeline has been replaced with MetaWorks: A flexible, scalable bioinformatic pipeline for multi-marker biodiversity assessments available from https://github.com/terrimporter/MetaWorks**

This repository outlines how 16S rDNA metabarcodes are processed by Teresita M. Porter. **SCVUR** refers to the programs, algorithms, and reference datasets used in this data flow: **S**EQPREP, **C**UTADAPT, **V**SEARCH, **U**SEARCH-UNOISE, **R**DP Classifier. 

## Overview

This pipeline begins with raw paired-end Illumina MiSeq fastq.gz files.  Reads are paired.  Primers are trimmed.  All the samples are pooled for a global analysis.  Reads are dereplicated and denoised producing a reference set of exact sequence variants (ESVs).  These ESVs are taxonomically assigned using the 16S reference set available with the RDP Classifier (Wang et al., 2007) available from https://sourceforge.net/projects/rdp-classifier/ .

This data flow has been developed using a conda environment and snakemake pipeline for improved reproducibility.  It will be updated on a regular basis so check for the latest version at https://github.com/terrimporter/SCVUR_16S_metabarcode_pipeline/releases

## Outline

[How to cite](#How-to-cite)  

[Standard pipeline](#standard-pipeline) 

[Implementation notes](#implementation-notes)  

[References](#references)  

[Acknowledgements](#acknowledgements)  

## How to cite

You can cite this repository directly:  
Teresita M. Porter. (2019, January 23). SCVUR 16S Metabarcode Pipeline (Version v3.0.2). Zenodo. http://doi.org/10.5281/zenodo.4741499  

## Pipeline details

If you are comfortable reading code, read through the snakefile to see how the pipeline runs as well as which programs and versions are used. Otherwise you can just list all the programs in the conda environment, see [Implementation notes](#implementation-notes). 

Raw paired-end reads are merged using SEQPREP v1.3.2 from bioconda (St. John, 2016).  This step looks for a minimum Phred quality score of 20 in the overlap region, requires at least 25bp overlap.

Primers are trimmed in two steps using CUTADAPT v2.6 from bioconda (Martin, 2011).  This step looks for a minimum Phred quality score of at least 20 at the ends, forward primer is trimmed first, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.  The output from the first step, is used as in put for the second step.  This step looks for a minimum Phred quality score of at least 20 at the ends, the reverse primer is trimmed, no more than 3 N's allowed, trimmed reads need to be at least 150 bp, untrimmed reads are discarded.

Files are reformatted and samples are combined for a global analysis.

Reads are dereplicated (only unique sequences are retained) using VSEARCH v2.13.6 from bioconda (Rognes et al., 2016).

Denoised exact sequence variants (ESVs) are generated using VSEARCH with the unoise3 algorithm (Edgar, 2016).  This step removes any PhiX contamination, sequences with predicted errors, and rare sequences.  This step also produces zero-radius OTUs (Zotus) also referred to commonly as amplicon sequence variants (ASVs), ESVs, or 100% operational taxonomic unit (OTU) clusters.  Here, we define rare sequences to be sequence clusters containing only one or two reads (singletons and doubletons) and these are removed as 'noise'.  Putative chimeric sequences are then removed using the uchime3_denovo algorithm in VSEARCH.

Primer trimmed reads (cat.fasta2.gz) are mapped to the denoised ESVs (cat.denoised.nonchimeras) in VSEARCH using the --search_exact method to generate the ESV.table .

16S taxonomic assignments are made using the Ribosomal Database classifier v2.12 (RDP classifier) available from https://sourceforge.net/projects/rdp-classifier/ (Wang et al., 2007).

The final output is reformatted to add read numbers from the ESV.table for each sample and column headers to improve readability (rdp.csv).  This final output file, rdp.csv, contains all ESVs for each sample, associated read counts, as well as taxonomic assignments with bootstrap support values.  rdp.csv can be read into R where the table can be filtered, reformatted, then reshaped to create a matrix of ESVs x samples containing read numbers to conduct standard biodiversity analyses.

Read and ESV statistics are provided for various steps of the program are also provided.

### Prepare your environment to run the pipeline

1. This pipeline includes a conda environment that provides most of the programs needed to run this pipeline (SNAKEMAKE, SEQPREP, CUTADAPT, VSEARCH, etc.).

```linux
# Create the environment from the provided environment.yml file
conda env create -f environment.yml

# Activate the environment
conda activate myenv.3
```

2. The pipeline also requires the RDP classifier for the taxonomic assignment step.  Although the RDP classifier v2.2 is available through conda, a newer v2.12 is available form SourceForge at https://sourceforge.net/projects/rdp-classifier/ .  Download it and take note of where the classifier.jar file is as this needs to be added to config.yaml .

The RDP classifier comes with the training sets to classify 16S (default), fungal ITS and fungal LSU rDNA sequences.

```linux
RDP:
    jar: "/path/to/rdp_classifier_2.12/dist/classifier.jar"
	c: 0
	f: "fixrank"
```

3. In most cases, your raw paired-end Illumina reads can go into a directory called 'data' which should be placed in the same directory as the other files that come with this pipeline.

```linux
# Create a new directory to hold your raw data
mkdir data
```

4. Please go through the config.yaml file and edit directory names, filename patterns, etc. as necessary to work with your filenames.

5. Be sure to edit the first line of each Perl script (shebang) in the perl_scripts directory to point to where Perl is installed.

```linux
# The usual shebang if you already have Perl installed
#!/usr/bin/perl

# Alternate shebang if you want to run perl using the conda environment (edit this)
#!/path/to/miniconda3/envs/myenv.3/bin/perl
```

### Run the pipeline

Run snakemake by indicating the number of jobs or cores that are available to run the whole pipeline.  

```linux
snakemake --jobs 24 --snakefile snakefile --configfile config.yaml
```

When you are done, deactivate the conda environment:

```linux
conda deactivate
```

## Implementation notes

### Installing Conda and Snakemake

Conda is an open source package and envirobnment management system.  Miniconda is a lightweight version of conda that only contains conda, python, and their dependencies.  Using conda and the environment.yml file provided here can help get all the necessary programs in one place to run this pipeline.  Snakemake is a Python-based workflow management tool meant to define the rules for running this bioinformatic pipeline.  There is no need to edit the snakefile or snakefile_alt files directly.  Changes to select parameters can be made in the config.yaml pipeline.  If you install Conda and activate the environment provided, then you will also get the correct versions of the open source programs used in this pipeline including Snakemake v3.13.3.

Install miniconda as follows:

```linux
# Download miniconda3
wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh

# Install miniconda3
sh Miniconda3-latest-Linux-x86_64.sh

# Add conda to your PATH, ex. to ~/bin
cd ~/bin
ln -s miniconda3/bin/conda conda
```

### Check program versions

Ensure that the correct programs from the environment are being used.

```linux
# create conda environment from file
conda env create -f environment.yml

# activate the environment
conda activate myenv

# list all programs available in the environment at once
conda list > programs.list

# or, inidivdually check that key programs in the conda environment are being used
which SeqPrep
which cutadapt
which vsearch
which perl

# then, check their version numbers one at a time
cutadapt --version
vsearch --version
```

Version numbers are also tracked in the snakefile.

Note that commercial software (ex. USEARCH) and programs not available from conda need to be installed on your system and executable in your PATH (see [Standard pipeline](#standard-pipeline) "Prepare your environment to run the pipeline").

### Batch renaming of files

Sometimes it is necessary to rename large numbers of sequence files.  I prefer to use Perl-rename (Gergely, 2018) that is available at https://github.com/subogero/rename as opposed to linux rename.  I prefer the Perl implementation so that you can easily use regular expressions.  I first run the command with the -n flag so you can review the changes without making any actual changes.  If you're happy with the results, re-run without the -n flag.

```linux
rename -n 's/PATTERN/NEW PATTERN/g' *.gz
```

### Symbolic links

Instead of continually traversing nested directories to get to files, I create symbolic links to target directories in a top level directory.  Symbolic links can also be placed in your ~/bin directory that point to scripts that reside elsewhere on your system.  So long as those scripts are executable (e.x. chmod 755 script.plx) then the shortcut will also be executable without having to type out the complete path or copy and pasting the script into the current directory.  This can be especially useful so that you don't have to maintain multiple copies of large raw read files in different places.

```linux
ln -s /path/to/target/directory shortcutName
ln -s /path/to/script/script.sh commandName
```

## References

Edgar, R. C. (2016). UNOISE2: improved error-correction for Illumina 16S and ITS amplicon sequencing. BioRxiv. doi:10.1101/081257 .  Available from: https://www.drive5.com/ 

Gergely, S. (2018, January). Perl-rename. Retrieved from https://github.com/subogero/rename  

Martin, M. (2011). Cutadapt removes adapter sequences from high-throughput sequencing reads. EMBnet. Journal, 17(1), pp–10.  Available from: http://cutadapt.readthedocs.io/en/stable/index.html

Pruesse E, Quast C, Knittel K, Fuchs BM, Ludwig WG, Peplies J, Glöckner FO (2007) SILVA: a comprehensive online resource for quality checked and aligned ribosomal RNA sequence data compatible with ARB. Nucl. Acids Res. 35:7188-7196

Rognes, T., Flouri, T., Nichols, B., Quince, C., & Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. PeerJ, 4, e2584. doi:10.7717/peerj.2584  

St. John, J. (2016, Downloaded). SeqPrep. Retrieved from https://github.com/jstjohn/SeqPrep/releases 

Tange, O. (2011). GNU Parallel - The Command-Line Power Tool. ;;Login: The USENIX Magazine, February, 42–47.  Available from: https://www.gnu.org/software/parallel/

Wang, Q., Garrity, G. M., Tiedje, J. M., & Cole, J. R. (2007). Naive Bayesian Classifier for Rapid Assignment of rRNA Sequences into the New Bacterial Taxonomy. Applied and Environmental Microbiology, 73(16), 5261–5267. doi:10.1128/AEM.00062-07 .  Available from: https://sourceforge.net/projects/rdp-classifier/

## Acknowledgements

I would like to acknowedge funding from the Canadian government through the Genomics Research and Development Initiative (GRDI) EcoBiomics project.

Last updated: May 6, 2021
