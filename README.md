Introduction
------------
transIndel is used to detect indels (insertions and deletions) from DNA-seq or RNA-seq data by parsing chimiric alignments from BWA-MEM. 

Prerequisites
----------------
Samtools/0.1.18 (http://samtools.sourceforge.net/)
Python packages:
* Pysam/0.7.7 or newer (https://code.google.com/p/pysam/)
* HTSeq/0.6.1 or newer (https://pypi.python.org/pypi/HTSeq)

Getting Soure Code
------------------
	git clone git://github.com/cauyrd/transIndel.git
	cd ScanIndel
Running transIndel 
-----------------
#### STEP 1: Build new BAM file with redefined CIGAR string
* analyzing DNA-seq data (whole genome seq/exome-seq/targeted capture)
	```
	python transIndel_build_DNA.py -i input_bam_file -o output_bam_file [options]
	```
	
* analyzing RNA-seq data 
	```
	python transIndel_build_RNA.py -i input_bam_file -r reference_genome_fasta -g gtf_file -o output_bam_file [options]
	```
#### Options:
	```
	--mapq_cutoff				:minimal MapQ in SAM for supporting reads (default 15)
	--max_del_length			:maximum deletion length to be detected (default 1Mbp)
	-h --help					:produce this menu
	-v --version				:show version of this tool
	```
#### Input:
	```
	input_bam_file   						:input BAM file is produced by BWA-MEM and is sorted and indexed.
	reference_genome_fasta (for RNA-seq)    :reference genome in FastA format
	gtf_file (for RNA-seq)    				:gene annotation file in GTF format
	```
#### Output:
	```
	your_output_bam_file		:BAM file for CIGAR string redefinement.
	```
	transIndel generates the following optional fields in output BAMs



	Tag | Meaning
	------------ | -------------
	OA | original representative alignment; format: (pos,CIGAR)
	JM | splicing junction reads; infered from GTF or splicing motif (used in RNA-seq BAM)
	

#### STEP 2: Call indel
* Option 1: using transIndel_call.py script
	```
	python transIndel_call.py -i input_bam_from_transIndel_build -o output_vcf_filename_prefix [options]	
	```
#### Options:
	```
	 -c							:minimal observation count for Indel (default 4)
	 -d							:minimal depth to call Indel (default 10)
	 -f							:minimal variant allele frequency (default 0.1)
	 -l							:minimal indel length to report (default 10)
	 -m							:minimal mapq of read from BAM file to call indel (default 15)
	 -t							:Limit analysis to targets listed in the BED-format FILE or a samtools region string
	 -h --help					:produce this menu
	 -v --version				:show version of this tool
	 ```
#### Input:
	```
	input_bam_file   			:input BAM file is produced by transIndel_build.py
	```
#### Output:
	```
	output_vcf_file   			:Reported Indels with VCF format
	```
	
* Option 2: using existing variant caller (e.g. freebayes, GATK)
	```
	following the specific variant caller's manual
	```
