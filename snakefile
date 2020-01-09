""" SCVUR v3 16S metabarcode pipeline

Metabarcode pipeline to process Illumina MiSeq reads as follows:

1_Pair reads
2_Trim reads
3_Concatenate samples for global analysis
4_Dereplicate reads
5_Denoise ESVs
6_Remove chimeras
7_Create ESV x sample table
8_Taxonomic assignment

"""

#######################################################################
# Read in vars from config.yaml
# Generate new vars, dirs, and filenames
# All code not in rules (this section) is run first before any rules (next section onwards)

# 1_Pair reads
raw_sample_read_wildcards=config["raw_sample_read_wildcards"]
print(raw_sample_read_wildcards)

SAMPLES,READS = glob_wildcards(raw_sample_read_wildcards)
SAMPLES_UNIQUE = list(set(SAMPLES))

raw_sample_forward_wildcard=config["raw_sample_forward_wildcard"]

dir=config["dir"]
R1stats_out="{sample}.R1stats"
raw1_stats_out=dir+"/stats/"+R1stats_out
cat_raw1_stats_out=dir+"/stats/R1.stats"

raw_sample_reverse_wildcard=config["raw_sample_reverse_wildcard"]
R2stats_out="{sample}.R2stats"
raw2_stats_out=dir+"/stats/"+R2stats_out
cat_raw2_stats_out=dir+"/stats/R2.stats"

fastq_gz="{sample}.fastq.gz"
seqprep_out=dir+"/paired/"+fastq_gz
Pstats_out="{sample}.Pstats"
paired_stats_out=dir+"/stats/"+Pstats_out
cat_paired_stats_out=dir+"/stats/paired.stats"

# 2_Trim reads
cutadapt_f_out=dir+"/Fprimer_trimmed/"+fastq_gz
Fstats_out="{sample}.Fstats"
Ftrimmed_stats_out=dir+"/stats/"+Fstats_out
cat_Ftrimmed_stats_out=dir+"/stats/Ftrimmed.stats"

fasta_gz="{sample}.fasta.gz"
cutadapt_r_out=dir+"/Rprimer_trimmed/"+fasta_gz
Rstats_out="{sample}.Rstats"
Rtrimmed_stats_out=dir+"/stats/"+Rstats_out
cat_Rtrimmed_stats_out=dir+"/stats/Rtrimmed.stats"

# 3_Concatenate samples for global analysis
fasta="{sample}.fasta"
concatenate_pattern=dir+"/"+fasta
output=dir+"/cat.fasta"
output2=dir+"/cat.fasta2"
gzip_out=dir+"/cat.fasta2.gz"

# 4_Dereplicate reads
dereplicate_out=dir+"/cat.uniques"
dereplicate_log=dir+"/dereplication.log"

# 5_Create denoised ESVs
denoise_out=dir+"/cat.denoised"
denoise_log=dir+"/usearch.log"

# 6_Chimera check
chimera_out=dir+"/cat.denoised.nonchimeras"
chimera_log=dir+"/chimeraRemoval.log"

# 7_Create ESV table
table_out=dir+"/ESV.table"
table_log=dir+"/table.log"

# 8_Taxonomic assignment
rdp_out=dir+"/rdp.out.tmp"
rdp_csv=dir+"/rdp.csv.tmp"
rdp_csv2=dir+"/rdp.csv"

#######################################################################
# This rule defines the final target file that is needed from the pipeline
# Hashed lines (begin with # symbol) are not executed
# Target files need to be separated with a comma
# Snakemake looks at all input and output files to determine which rules need to be run to create the target file

rule all:
	input:
		# Rule testing [In order of execution]:
		# Calculate raw stats
#		expand(raw1_stats_out, sample=SAMPLES_UNIQUE),
#		expand(raw2_stats_out, sample=SAMPLES_UNIQUE),
		cat_raw1_stats_out,
		cat_raw2_stats_out,
		# 1_Pair reads
#		expand(seqprep_out, sample=SAMPLES)
		# Calculate paired stats
		expand(paired_stats_out, sample=SAMPLES_UNIQUE),
		cat_paired_stats_out,
		# 2_Trim reads (forward)
#		expand(cutadapt_f_out, sample=SAMPLES)
		# Calculate forward trimmed stats
		expand(Ftrimmed_stats_out, sample=SAMPLES_UNIQUE),
		cat_Ftrimmed_stats_out,
		# 2_Trim reads (reverse)
#		expand(cutadapt_r_out, sample=SAMPLES)
		# Calculate reverse trimmed stats
		expand(Rtrimmed_stats_out, sample=SAMPLES_UNIQUE),
		cat_Rtrimmed_stats_out,
		# 2_Trim reads (edit fasta header)
#		expand(concatenate_pattern, sample=SAMPLES_UNIQUE)
		# 3_Concatenate samples for global analysis
#		output
		# 3_Concatenate samples for global analysis (edit fasta header)
#		output2
		# 3_Concatenate samples for global analysis (compress file)
#		gzip_out
		# 4_Dereplicate reads
#		dereplicate_out
		# 5_Denoise reads
#		denoise_out
		# 6_Chimera check
#		chimera_out
		# 7_Create ESV table
#		table_out
		# 8_Taxonomic assignment
#		rdp_out
		# 8_Taxonomic assignment (add read numbers)
#		rdp_csv
		# 8_Taxonomic assignment (edit ESV id's to include amplicon name) [Final output file]
		rdp_csv2
	 
#######################################################################
# Count raw forward reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > R1.stats'

rule raw1_stats:
	input: 
		raw_sample_forward_wildcard
	output:
		temp(raw1_stats_out)
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate R1 reads

rule cat_raw1_stats:
	input:
		expand(raw1_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_raw1_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Count raw reverse reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > R2.stats'

rule raw2_stats:
	input: 
		raw_sample_reverse_wildcard
	output:
		temp(raw2_stats_out)
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate R2 reads

rule cat_raw2_stats:
	input:
		expand(raw2_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_raw2_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Pair forward and reverse reads with SeqPrep

rule pair_reads:
	version: "1.3.2"
	input:
		f=raw_sample_forward_wildcard,
		r=raw_sample_reverse_wildcard
	output:
		X1=temp("{sample}_R1.out"),
		X2=temp("{sample}_R2.out"),
		s=seqprep_out
	shell:
		"SeqPrep -f {input.f} -r {input.r} -1 {output.X1} -2 {output.X2} -q {config[SEQPREP][q]} -s {output.s} -o {config[SEQPREP][o]}"

#######################################################################
# Count paired reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > paired.stats'

rule paired_stats:
	input: 
		seqprep_out
	output:
		temp(paired_stats_out)
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate paired reads

rule cat_paired_stats:
	input:
		expand(paired_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_paired_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Trim forward primer with CUTADAPT

rule trim_forward_primer:
	version: "2.6"
	input:
		seqprep_out
	output:
		cutadapt_f_out
	shell:
		"cutadapt -g {config[CUTADAPT_FWD][g]} -m {config[CUTADAPT_FWD][m]} -q {config[CUTADAPT_FWD][q]} --max-n={config[CUTADAPT_FWD][mn]} --discard-untrimmed -o {output} {input}"

#######################################################################
# Count forward primer trimmed reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > Ftrimmed.stats'

rule Ftrimmed_stats:
	input: 
		cutadapt_f_out
	output:
		temp(Ftrimmed_stats_out)
	shell:
		"perl perl_scripts/fastq_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate Ftrimmed reads

rule cat_Ftrimmed_stats:
	input:
		expand(Ftrimmed_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_Ftrimmed_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Trim reverse primer with CUTADAPT

rule trim_reverse_primer:
	version: "2.6"
	input:
		cutadapt_f_out
	output:
		cutadapt_r_out
	shell:
		"cutadapt -a {config[CUTADAPT_REV][a]} -m {config[CUTADAPT_REV][m]} -q {config[CUTADAPT_REV][q]} --max-n={config[CUTADAPT_REV][mn]} --discard-untrimmed -o {output} {input}"

#######################################################################
# Count reverse primer trimmed reads
# For each file calculates number of reads (total), length of reads (min, max, mean, median, mode)
# To summarize output, go to directory and enter 'cat *.stats' > Rtrimmed.stats'

rule Rtrimmed_stats:
	input: 
		cutadapt_r_out
	output:
		temp(Rtrimmed_stats_out)
	shell:
		"perl perl_scripts/fasta_gz_stats.plx {input} >> {output}"

#######################################################################
# Concatenate Rtrimmed reads

rule cat_Rtrimmed_stats:
	input:
		expand(Rtrimmed_stats_out, sample=SAMPLES_UNIQUE)
	output:
		cat_Rtrimmed_stats_out
	shell:
		"cat {input} >> {output}"

#######################################################################
# Edit fasta header with Perl script

rule edit_fasta_header1:
	input:
		cutadapt_r_out
	output:
		temp(concatenate_pattern)
	threads: 1
	shell:
		"perl perl_scripts/rename_fasta_gzip.plx {input} > {output}"

#######################################################################
# Concatenate all samples into single file for global analysis using Linux cat
# Only use a single thread to work right

rule concatenate:
	input:
		expand(concatenate_pattern, sample=SAMPLES_UNIQUE)
	output:
		temp(output)
	threads: 1
	shell:
		"cat {input} >> {output}"

#######################################################################
# Edit fasta header again using GNU sed

rule edit_fasta_header2:
	version: "4.7"
	input:
		output
	output:
		temp(output2)
	shell:
		"sed -e 's/-/_/g' {input} >> {output}"

#######################################################################
# Compress into smaller file using Linux gzip

rule compress:
	version: "1.7"
	input:
		output2
	output:
		gzip_out
	shell:
		"gzip -c {input} > {output}"

#######################################################################
# Dereplicate and track reads with VSEARCH

rule dereplicate:
	version: "2.14.1"
	input:
		gzip_out
	output:
		dereplicate_out
	log: dereplicate_log
	shell:
		"vsearch --derep_fulllength {input} --output {output} --sizein --sizeout --log {log}"
	
#######################################################################
# Denoise with VSEARCH using unoise3 algorithm (chimera removal must be done separately in next step)

rule denoise:
	version: "2.14.1"
	input:
		dereplicate_out
	output:
		denoise_out
	log: denoise_log
	shell:
		"vsearch --cluster_unoise {input} --sizein --sizeout --minsize {config[VSEARCH_DENOISE][minsize]} --centroids {output} --log {log}"

#######################################################################
# Chimera removal with VSEARCH using uchime_denovo3 algorithm

rule chimera_removal:
	version: "2.14.1"
	input:
		denoise_out
	output:
		chimera_out
	log: chimera_log
	shell:
		"vsearch --uchime3_denovo {input} --sizein --xsize --nonchimeras {output} --relabel 'Zotu' --log {log}"
	
#######################################################################
# Create ESV table

rule create_ESV_table:
	version: "2.14.1"
	input:
		vsearch_global=gzip_out,
		db=chimera_out
	output:
		table_out
	threads: config["VSEARCH_TABLE"]["t"]
	log: table_log
	shell:
		"vsearch --threads {threads} --usearch_global {input.vsearch_global} --db {input.db} --id 1.0 --otutabout {output} --log {log}"

#######################################################################
# Taxonomic assignment

# Use the RDP classifier for taxonomic assignment
# Make sure this is installed locally (do not use the old 2.2 version available from conda), download the v2.12 from sourceforge, indicate the path to the classifeir in the config file
# Be sure to indicate the location of the trained reference database properties file in the config file

rule taxonomic_assignment:
	version: "2.12"
	input:
		chimera_out
	output:
		temp(rdp_out)
	shell:
		"java -Xmx8g -jar {config[RDP][jar]} classify -c {config[RDP][c]} -f {config[RDP][f]} -o {output} {input}"
	
#######################################################################
# Map read number

# Edit path in top line of Perl script to match path to perl in conda environment

rule map_read_number:
	input:
		table=table_out,
		rdp=rdp_out
	output:
		temp(rdp_csv)
	shell:
		"perl perl_scripts/add_abundance_to_rdp_out4b.plx {input.table} {input.rdp} > {output}"

#######################################################################
# Edit csv to add amplicon name to Zotu
# This is needed when data from different amplicons are combined, so that the Zotu's remain unique

rule edit_csv:
	input:
		rdp_csv
	output:
		rdp_csv2
	shell:
		"sed -e {config[SED]} {input} > {output}"


# End of Metabarcode pipeline
