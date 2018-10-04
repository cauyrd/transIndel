#!/usr/bin/python
#-*- coding: utf-8 -*-
#===============================================================================
#
#         FILE: transIndel_build_DNA.py
#
#        USAGE: transIndel_build_DNA.py -h
#
#  DESCRIPTION: correct CIGAR string for indel
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
try: import numpy as np
except: sys.exit('numpy module not found.\nPlease install it before.')

__version__ = 'v1.0'

def detect_softclip_mode(cigarstring):
	# if cigarstring is 40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')] with the statement below
	# mode 1 = MS; 2=SM; 3=other
	if 'I' in cigarstring or 'D' in cigarstring or 'H' in cigarstring:
		return 3,-1,-1
	cigartuple = re.findall(r'(\d+)(\w)', cigarstring)
	if cigartuple[0][1] == 'S':
		if cigartuple[-1][1] == 'S':
			if int(cigartuple[0][0]) > int(cigartuple[-1][0]):
				return 2,int(cigartuple[1][0]),int(cigartuple[-1][0])
			else:
				return 1,int(cigartuple[1][0]),int(cigartuple[0][0])
		else:
			return 2,int(cigartuple[1][0]),0
	else:
		return 1,int(cigartuple[0][0]),0

def detect_sv_from_cigar(chr, read, mapq_cutoff, max_del_len):
	strand = '-' if read.is_reverse else '+'
	try:
		chimeric_aln = read.get_tag('SA').split(';')[:-1]
	except KeyError:
		return 'NA', [] 
	if len(chimeric_aln) == 1:
		chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.get_tag('SA').split(',')
		if int(mapq_sa) < mapq_cutoff:
			return 'NA', []
	else:
		return 'NA', []
	read_mode,read_match,read_clipped = detect_softclip_mode(read.cigarstring)
	sa_mode,sa_match,sa_clipped = detect_softclip_mode(cigar_sa)
	if read_mode == 3 or sa_mode == 3:
		return 'NA', []
	target_start = 0
	target_end = 0
	newcigar = []
	newpos = -1
	if chr == chr_sa:
		if strand_sa == strand: # deletion, insertion, duplication
			if read_mode == 2 and sa_mode == 1:
				target_start = int(pos_sa) - 1
				target_end = read.aend
				if sa_clipped != 0:
					newcigar.append([4,sa_clipped])
				newcigar.append([0,sa_match])
				target_offset = target_end - target_start
				query_offset = read.rlen - sa_clipped - read_clipped
				indel_size = query_offset - target_offset
				newpos = target_start
				if indel_size == 0: 
					return 'NA', []
				elif indel_size < 0: # deletion
					if abs(indel_size) > max_del_len:
						return 'NA', []
					newcigar.append([2,-1*indel_size])
					newcigar.append([0,query_offset-sa_match])
				elif indel_size >= query_offset: 
					return 'NA', []
				else: # insertion
					newcigar.append([1,indel_size])
					last_match_size = read.aend - (int(pos_sa)-1+sa_match)
					if last_match_size <= 0:
						for each in newcigar:
							if each[0] == 0:
								each[1] = each[1]+last_match_size
								break
					else:
						newcigar.append([0,last_match_size])
				if read_clipped != 0:
					newcigar.append([4,read_clipped])
			elif read_mode == 1 and sa_mode == 2:
				target_start = read.pos
				target_end = int(pos_sa) -1 + sa_match
				if read_clipped != 0:
					newcigar.append([4,read_clipped])
				newcigar.append([0,read_match])
				target_offset = target_end - target_start
				query_offset = read.rlen - sa_clipped - read_clipped
				indel_size = query_offset - target_offset
				newpos = target_start
				if indel_size == 0:
					return 'NA', []
				elif indel_size < 0:
					if abs(indel_size) > max_del_len:
						return 'NA', []
					newcigar.append([2,-1*indel_size])
					newcigar.append([0,query_offset-read_match])
				elif indel_size >= query_offset: 
					return 'NA', []
				else:
					newcigar.append([1,indel_size])
					last_match_size = int(pos_sa)-1+sa_match -  read.aend
					if last_match_size <= 0:
						for each in newcigar:
							if each[0] == 0:
								each[1] = each[1]+last_match_size
								break
					else:
						newcigar.append([0,last_match_size])
				if sa_clipped != 0:
					newcigar.append([4,sa_clipped])
			else:
				return 'NA', []
		else:
			return 'NA', []
	else:
		return 'NA', []
	return map(tuple,newcigar), newpos

def softclipping_realignment(mapq_cutoff, max_del_len, input, output):
	bwa_bam = pysam.Samfile(input,'rb')
	output_bam = pysam.Samfile(output+'.temp.bam','wb', template=bwa_bam)
	try:
		for read in bwa_bam.fetch(until_eof=True):
			if read.mapq >= mapq_cutoff and not read.is_secondary and not read.has_tag('XA'): 
				chr = bwa_bam.getrname(read.rname)
				newcigar,newpos =  detect_sv_from_cigar(chr,read,mapq_cutoff,max_del_len)
				if newcigar != 'NA' and newcigar != read.cigar: 
					read.setTag('OA', str(read.pos+1)+','+read.cigarstring)
					read.cigar, read.pos = newcigar, newpos
			output_bam.write(read)
	except ValueError as e:
		print >> sys.stderr, 'Bam index file is not found!', e
		sys.exit(1)
	bwa_bam.close()
	output_bam.close()
	try:
		subprocess.check_call("samtools sort {0}.temp.bam -o {0}".format(output), shell=True)
	except subprocess.CalledProcessError as e:
		print >> sys.stderr, 'Execution failed for samtools:', e
		sys.exit(1)
	subprocess.check_call("samtools index {}".format(output), shell=True)

def usage():
	"""helping information"""
	print 'Usage:'
	print ' python transIndel_build_DNA.py -i input_bam_file -o output_bam_file[opts]'
	print 'Opts:'
	print ' --mapq_cutoff  :minimal MapQ in SAM for support SV event, default 15'
	print ' --max_del_length  :maximum deletion length to be detected (10e6)'
	print ' -h --help :produce this menu'
	print ' -v --version :show version of this tool'
	print 'author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'
	print 'version: ' + __version__

def external_tool_checking():
	"""checking dependencies are installed"""
	software = ['samtools']
	cmd = "which"
	for each in software:
		try:
			path = subprocess.check_output([cmd, each], stderr=subprocess.STDOUT)
		except subprocess.CalledProcessError:
			print >> sys.stderr, "Checking for '" + each + "': ERROR - could not find '" + each + "'"
			print >> sys.stderr, "Exiting."
			sys.exit(0)
		print "Checking for '" + each + "': found " + path
	return path.rstrip()

def version_checking(samtools_path):
	'''HTSlib-based Samtools is needed!
	'''
	try:
		response = subprocess.check_output('{} --version-only'.format(samtools_path), shell=True, stderr=subprocess.STDOUT)
	except subprocess.CalledProcessError:
		print >> sys.stderr, "HTSlib-based Samtools is not available, please install latest Samtools at http://www.htslib.org/"
		print >> sys.stderr, "Exiting."
		sys.exit(0)


def main():

	# parameters parsing
	input = ''
	output = ''
	mapq_cutoff = 15 
	max_del_len = 10e6
	try:
		opts, args = getopt.getopt(sys.argv[1:], 'i:o:h:v', ['mapq_cutoff=', 'max_del_length=', 'help', 'version'])
		if not opts:
			print "Please use the -h or --help option to get usage information."
			sys.exit(0)
	except getopt.GetoptError as err:
		print >> sys.stderr, err
		usage()
		sys.exit(2)
	for o, a in opts:
		if o == '-i': input = a
		elif o == '-o': output = a
		elif o == '--mapq_cutoff': mapq_cutoff = int(a)
		elif o == '--max_del_length': max_del_len = int(a)
		elif o in ('-v', '--version'):
			print 'transIndel ' + __version__
			sys.exit(0)
		elif o in ('-h', '--help'):
			usage()
			sys.exit(0)
		else:
			assert False, "unhandled option"
	if not input or not output:
		print >> sys.stderr, 'Please specify -i and -o for input and output files.'
		usage()
		sys.exit(1)

	print 'transIndel build starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S")
	start = time.time()

	#check external tools used
	path = external_tool_checking()
	version_checking(path)
	# CIGAR string refinement or add SV tag 
	softclipping_realignment(mapq_cutoff, max_del_len, input, output)

	os.system('rm '+output+'.temp.*')
	print "transIndel build running done: " + time.strftime("%Y-%m-%d %H:%M:%S")
	end = time.time()
	print 'transIndel build takes ' + str(end - start) + ' seconds.'

if __name__ == '__main__':
	main()
