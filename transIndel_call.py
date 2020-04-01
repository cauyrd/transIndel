#!/usr/bin/env python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: transIndel_call.py
#
#        USAGE: transIndel_call.py -h
#
#  DESCRIPTION: call indel from transIndel_build.py output
#
#      OPTIONS: ---
# REQUIREMENTS: See README.md file
#         BUGS: ---
#        NOTES: ---
#       AUTHOR: Rendong Yang (yang4414@umn.edu)
# ORGANIZATION:
#===============================================================================
import sys, os, time, subprocess, getopt, re
try: import pysam
except: sys.exit('pysam module not found.\nPlease install it before.')
#try: import numpy as np
#except: sys.exit('numpy module not found.\nPlease install it before.')

__version__ = 'v0.1'

def vcf_header(output_prefix):
	header = ['##fileformat=VCFv4.1']
	header.append('##source=transIndel '+__version__)
	header.append('##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
	header.append('##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">')
	header.append('##INFO=<ID=AO,Number=1,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">')
	header.append('##INFO=<ID=AB,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads">')
	header.append('##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event, INS, DEL">')
	header.append('##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
	header.append('##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation (under development)">')
	header.append('##INFO=<ID=END,Number=1,Type=Integer,Description="end position of the indel">')
	header.append('##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect indel (under development)">')
	header.append('##ALT=<ID=DEL,Description="Deletion">')
	header.append('##ALT=<ID=INS,Description="Insertion">')
	header.append('##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
	header.append('#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t'+output_prefix)
	return '\n'.join(header)

def is_sv(read,ref_site,mapq_cutoff):
	if read.alignment.has_tag('SV'):
		chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.alignment.get_tag('SA').split(',')
		if int(mapq_sa) < mapq_cutoff:
			return False
		sv_type, sv_pos, sv_length, read_mode, sa_mode = read.alignment.get_tag('SV').split(',') 
		if int(sv_pos) - 1 + (int(read_mode) - 2) == ref_site:
			return True
		elif sv_type != 'TRA' and int(sv_pos) - 1 + int(sv_length) + (int(read_mode) - 2) == ref_site:
			return True
	return False

def sv_scan(input_bam,output_prefix,ref_genome,ao_cutoff,dp_cutoff,vaf_cutoff,indel_len_cutoff,target,mapq_cutoff):
	sv_set = set()
	samfile = pysam.Samfile(input_bam,'rb')
	fastafile = pysam.Fastafile(ref_genome) if ref_genome else None 
	vcffile = open(output_prefix+'.indel.vcf','w',0)
	print >> vcffile, vcf_header(output_prefix)
	vcf_field_gt = 'GT\t1/1'	
	regions = []
	try:
		ifp = open(target)
		for line in ifp:
			items = line.rstrip().split()
			regions.append(items[0]+':'+items[1]+'-'+items[2])
	except IOError:
		if target == '':
			regions = [None]
		else:
			regions = [target]
	for each_region in regions:
		try:
			for col in samfile.pileup(region=each_region):
				dp = col.n
				sv_dict = {}
				for read in col.pileups:
					if read.alignment.mapq >= mapq_cutoff: 
						if read.indel >= indel_len_cutoff:
							sv_id = 'INS,'+str(col.pos+1)+','+str(read.indel)
							try:
								sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
							except KeyError:
								ref_seq = read.alignment.query_sequence[read.query_position] 
								alt_seq = read.alignment.query_sequence[read.query_position:read.query_position+read.indel+1] 
								sv_dict[sv_id] = [ref_seq, alt_seq, 1]
						elif -1*read.indel >= indel_len_cutoff:
							sv_id = 'DEL,'+str(col.pos+1)+','+str(-1*read.indel)
							try:
								sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
							except KeyError:
								ref_seq = fastafile.fetch(col.reference_name,col.pos+1,col.pos+2) if fastafile else '.'
								sv_dict[sv_id] = [ref_seq, '<DEL>', 1]
						elif is_sv(read, col.pos, mapq_cutoff):
							sv_id = ','.join(read.alignment.get_tag('SV').split(',')[:3])
							try:
								sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
							except KeyError:
								sv_dict[sv_id] = ['.', '<'+sv_id.split(',')[0]+'>', 1]
				for sv_id in sv_dict: # only consider one event per genomic position
					if col.reference_name+','+sv_id in sv_set:
						continue
					ao = sv_dict[sv_id][-1]
					vaf = round(float(ao)/dp,2)
					if dp >= dp_cutoff and ao >= ao_cutoff and vaf >= vaf_cutoff:
						sv_set.add(col.reference_name+','+sv_id)
						sv_type, sv_pos, sv_length = sv_id.split(',')
						if sv_type == 'INS':
							end = int(sv_pos) + 1
						else:
							end = int(sv_pos) + int(sv_length)
						chr2 = col.reference_name if sv_type != 'TRA' else sv_length.split(':')[0]
						sv_length = sv_length if sv_type != 'TRA' else '0'
						vcf_field_info = 'NS=1;AO='+str(ao)+';DP='+str(dp)+';AB='+str(vaf)+';SVLEN='+sv_length+';SVTYPE='+sv_type+';SVMETHOD=transIndel_ALN;CHR2='+chr2+';END='+str(end)
						print >> vcffile, col.reference_name+'\t'+sv_pos+'\t.\t'+sv_dict[sv_id][0]+'\t'+sv_dict[sv_id][1]+'\t.\t.\t'+vcf_field_info+'\t'+vcf_field_gt
		except ValueError:
			continue
	samfile.close()
	return

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python transIndel_call.py -i input_bam_from_transIndel_build -o output_vcf_filename_prefix [opts]'
	print 'Opts:'
	print ' -r  :reference genome used for VCF REF column (required for valid VCF)'
	print ' -c  :minimal observation count for Indel (default 4)'
	print ' -d  :minimal depth to call Indel (default 10)'
	print ' -f  :minimal variant allele frequency (default 0.1)'
	print ' -l  :minimal indel length to report (default 10)'
	print ' -m  :minimal mapq of read from BAM file to call indel (default 15)'
	print ' -t  :Limit analysis to targets listed in the BED-format FILE or a samtools region string'
	print ' -h --help :produce this menu'
	print ' -v --version :show version of this tool'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: ' + __version__

def main():

	# parameters parsing
	input_bam = ''
	output_vcf_prefix = ''
	ref_genome = ''
	ao_cutoff = 4
	dp_cutoff = 10
	vaf_cutoff = 0.1
	indel_len_cutoff = 10 
	mapq_cutoff = 15 
	target = ''
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:o:r:c:d:f:l:m:t:h:v', ['help', 'version'])
		if not opts:
			print "Please use the -h or --help option to get usage information."
			sys.exit(0)
	except getopt.GetoptError as err:
		print >> sys.stderr, err
		usage()
		sys.exit(2)
	for o, a in opts:
		if o == '-i': input_bam = a
		elif o == '-o': output_vcf_prefix = a
		elif o == '-r': ref_genome = a
		elif o == '-c': ao_cutoff = int(a)
		elif o == '-d': dp_cutoff = int(a)
		elif o == '-f': vaf_cutoff = float(a)
		elif o == '-l': indel_len_cutoff = int(a)
		elif o == '-m': mapq_cutoff = int(a)
		elif o == '-t': target = a
		elif o in ('-v', '--version'):
			print 'transIndel' + __version__
			sys.exit(0)
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"
	if not input_bam or not output_vcf_prefix:
		print >> sys.stderr, 'Please specify input bam file and output vcf filename prefix.'
		usage()
		sys.exit(1)
	if indel_len_cutoff == 0:
		print >> sys.stderr, 'Minimal indel length has to be 1 base'
		usage()
		sys.exit(1)

	print 'transIndel call starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S")
	start = time.time()
	sv_scan(input_bam, output_vcf_prefix, ref_genome, ao_cutoff, dp_cutoff, vaf_cutoff, indel_len_cutoff, target, mapq_cutoff)
	print "transIndel call running done: " + time.strftime("%Y-%m-%d %H:%M:%S")
	end = time.time()
	print 'transIndel call takes ' + str(end - start) + ' seconds.'

if __name__ == '__main__':
	main()
