#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import time
import argparse
import re
import textwrap
try:
    import pysam
except:
    sys.exit('pysam module not found.\nPlease install it before.')

__version__ = 'v2.0'

def vcf_header(output_prefix):
    header = ['##fileformat=VCFv4.1']
    header.append('##source=transIndel '+__version__)
    header.append('##ALT=<ID=DEL,Description="Deletion">')
    header.append('##ALT=<ID=INS,Description="Insertion">')
    header.append(
        '##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">')
    header.append(
        '##INFO=<ID=DP,Number=1,Type=Integer,Description="Total read depth at the locus">')
    header.append(
        '##INFO=<ID=AO,Number=1,Type=Integer,Description="Alternate allele observations, with partial observations recorded fractionally">')
    header.append('##INFO=<ID=AB,Number=1,Type=Float,Description="Estimated allele frequency in the range (0,1], representing the ratio of reads showing the alternative allele to all reads">')
    header.append(
        '##INFO=<ID=SVTYPE,Number=1,Type=String,Description="The type of event, INS, DEL">')
    header.append(
        '##INFO=<ID=SVLEN,Number=1,Type=Integer,Description="Difference in length between REF and ALT alleles">')
    header.append(
        '##INFO=<ID=CHR2,Number=1,Type=String,Description="Chromosome for END coordinate in case of a translocation (under development)">')
    header.append(
        '##INFO=<ID=END,Number=1,Type=Integer,Description="end position of the indel">')
    header.append(
        '##INFO=<ID=SVMETHOD,Number=1,Type=String,Description="Type of approach used to detect indel (under development)">')
    header.append(
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">')
    header.append(
        '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t{}\n'.format(output_prefix))
    return '\n'.join(header)


def is_sv(read, ref_site, mapq_cutoff):
    if read.alignment.has_tag('SV'):
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.alignment.get_tag(
            'SA').split(',')
        if int(mapq_sa) < mapq_cutoff:
            return False
        sv_type, sv_pos, sv_length, read_mode, sa_mode = read.alignment.get_tag(
            'SV').split(',')
        if int(sv_pos) - 1 + (int(read_mode) - 2) == ref_site:
            return True
        elif sv_type != 'TRA' and int(sv_pos) - 1 + int(sv_length) + (int(read_mode) - 2) == ref_site:
            return True
    return False


def sv_scan(input_bam, output_prefix, ao_cutoff, dp_cutoff, vaf_cutoff, indel_len_cutoff, target, mapq_cutoff):
    sv_set = set()
    samfile = pysam.AlignmentFile(input_bam, 'rb')
    vcffile = open('{}.indel.vcf'.format(output_prefix), 'w', buffering=1)
    vcffile.write(vcf_header(output_prefix))
    vcf_field_gt = 'GT\t0/1'
    regions = []
    try:
        with open(target) as ifp:
            for line in ifp:
                items = line.rstrip().split()
                regions.append('{}:{}-{}'.format(items[0], items[1], items[2]))
    except TypeError:
        if target == '':
            regions = [None]
        else:
            regions = [target]
    for each_region in regions:
        try:
            for col in samfile.pileup(region=each_region):
                dp = col.nsegments
                sv_dict = {}
                for read in col.pileups:
                    if read.alignment.mapq >= mapq_cutoff:
                        if read.indel >= indel_len_cutoff:
                            sv_id = 'INS,'+str(col.pos+1)+','+str(read.indel)
                            try:
                                sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
                            except KeyError:
                                ref_seq = read.alignment.query_sequence[int(read.query_position)]
                                alt_seq = read.alignment.query_sequence[
                                    int(read.query_position):int(read.query_position)+read.indel+1]
                                sv_dict[sv_id] = [ref_seq, alt_seq, 1]
                        elif -1*read.indel >= indel_len_cutoff:
                            sv_id = 'DEL,'+str(col.pos+1) + \
                                ','+str(-1*read.indel)
                            try:
                                sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
                            except KeyError:
                                sv_dict[sv_id] = ['.', '<DEL>', 1]
                        elif is_sv(read, col.pos, mapq_cutoff):
                            sv_id = ','.join(
                                read.alignment.get_tag('SV').split(',')[:3])
                            try:
                                sv_dict[sv_id][-1] = sv_dict[sv_id][-1] + 1
                            except KeyError:
                                sv_dict[sv_id] = [
                                    '.', '<'+sv_id.split(',')[0]+'>', 1]
                for sv_id in sv_dict:  # only consider one event per genomic position
                    if col.reference_name+','+sv_id in sv_set:
                        continue
                    ao = sv_dict[sv_id][-1]
                    vaf = round(float(ao)/dp, 2)
                    if dp >= dp_cutoff and ao >= ao_cutoff and vaf >= vaf_cutoff:
                        sv_set.add(col.reference_name+','+sv_id)
                        sv_type, sv_pos, sv_length = sv_id.split(',')
                        if sv_type == 'INS':
                            end = int(sv_pos) + 1
                        else:
                            end = int(sv_pos) + int(sv_length)
                        chr2 = col.reference_name if sv_type != 'TRA' else sv_length.split(':')[
                            0]
                        sv_length = sv_length if sv_type != 'TRA' else '0'
                        vcf_field_info = 'NS=1;AO='+str(ao)+';DP='+str(dp)+';AB='+str(
                            vaf)+';SVLEN='+sv_length+';SVTYPE='+sv_type+';SVMETHOD=transIndel;CHR2='+chr2+';END='+str(end)
                        print(col.reference_name+'\t'+sv_pos+'\t.\t' + \
                            sv_dict[sv_id][0]+'\t'+sv_dict[sv_id][1] + \
                            '\t.\t.\t'+vcf_field_info+'\t'+vcf_field_gt, file=vcffile)
        except ValueError:
            sys.stderr.write('Error! {} at {}\n'.format(e, each_region))
            continue
    samfile.close()
    return True

def parse_args():
    class LengthAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if values <= 0:
                parser.error("Minimum length for {0} is 1".format(option_string))
            setattr(namespace, self.dest, values)

    parser = argparse.ArgumentParser(
        description = "%(prog)s -i input_bam_from_transIndel_build -o output_vcf_filename_prefix [opts]",
        epilog=textwrap.dedent('''transIndel:a splice-aware algorithm that detects the mid-sized insertions
            and large deletions using chimeric alignments from DNA-seq or RNA-seq data.
            Author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'''))
    parser.add_argument('-i', '--input', action='store', dest='input', help="Input BAM file", required=True)
    parser.add_argument('-o', '--output', action='store', dest='output', help="output VCF file prefix", required=True)
    parser.add_argument('-c', '--ao', action='store', dest='ao', type=int, help="minimal observation count for Indel (default: %(default)s)", default=4)
    parser.add_argument('-d', '--depth', action='store', dest='depth', type=int, help="minimal depth to call Indel (default: %(default)s)", default=10)
    parser.add_argument('-f', '--vaf', action='store', dest='vaf', type=float, help="minimal variant allele frequency (default: %(default)s)", default=0.1)
    parser.add_argument('-l', '--length', action=LengthAction, dest='length', type=int, help="minimal Indel length (>=1) to report (default: %(default)s)", default=10)
    parser.add_argument('-m', '--mapq', action='store', dest='mapq', type=int, help="minimal MAPQ of read from BAM file to call Indel (default: %(default)s)", default=15)
    parser.add_argument('-t', action='store', dest='region', help="Limit analysis to targets listed in the BED-format FILE or a samtools region string")
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))

    return parser


def main():
    if sys.version_info < (3,6):
        sys.exit('Sorry, this code need Python 3.6 or higher. Please update. Aborting...')
    parser = parse_args()
    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        options = parser.parse_args()

    print('transIndel call starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S"))
    start = time.time()
    sv_scan(options.input, options.output, options.ao, options.depth,
            options.vaf, options.length, options.region, options.mapq)
    print("transIndel call running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
    end = time.time()
    print('transIndel call takes ' + str(end - start) + ' seconds.')

if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
