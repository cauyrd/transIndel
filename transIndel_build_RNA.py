#!/usr/bin/env python
# -*- coding: utf-8 -*-
import sys
import os
import time
import subprocess
import argparse
import re
import textwrap
try:
    import pysam
except:
    sys.exit('pysam module not found.\nPlease install it before.')
try:
    import HTSeq
except:
    sys.exit('HTSeq module not found.\nPlease install it before.')

__version__ = 'v2.0'


def remove(infile):
    if os.path.isfile(infile):
        os.remove(infile)

def extract_splice_sites(file, bin):
    gtf_file = HTSeq.GFF_Reader(file)
    cvg = HTSeq.GenomicArray("auto", stranded=False, typecode='i')
    print("Loading GTF file...")
    for feature in gtf_file:
        if feature.type == 'exon':
            iv1 = HTSeq.GenomicInterval(
                feature.iv.chrom, feature.iv.start-bin, feature.iv.start+bin, '.')
            iv2 = HTSeq.GenomicInterval(
                feature.iv.chrom, feature.iv.end-bin, feature.iv.end+bin, '.')
            try:
                cvg[iv1] = 1
                cvg[iv2] = 1
            except IndexError:
                continue
    print("GTF file loaded!")
    return cvg


def detect_softclip_mode(cigarstring):
    # if cigarstring is 40M25N5M then cigartuple is [('40', 'M'), ('25', 'N'), ('5', 'M')] with the statement below
    # mode 1 = MS; 2=SM; 3=other
    if 'I' in cigarstring or 'D' in cigarstring or 'H' in cigarstring:
        return 3, -1, -1
    cigartuple = re.findall(r'(\d+)(\w)', cigarstring)
    if cigartuple[0][1] == 'S':
        if cigartuple[-1][1] == 'S':
            if int(cigartuple[0][0]) > int(cigartuple[-1][0]):
                return 2, int(cigartuple[1][0]), int(cigartuple[-1][0])
            else:
                return 1, int(cigartuple[1][0]), int(cigartuple[0][0])
        else:
            return 2, int(cigartuple[1][0]), 0
    else:
        if cigartuple[-1][1] == 'S':
            return 1, int(cigartuple[0][0]), 0
        else:
            return 3, -1, -1


def detect_sv_from_cigar(chrm, read, mapq_cutoff, max_del_len):
    strand = '-' if read.is_reverse else '+'
    try:
        chimeric_aln = read.get_tag('SA').split(';')[:-1]
    except KeyError:
        return 'NA', []
    if len(chimeric_aln) == 1:
        chr_sa, pos_sa, strand_sa, cigar_sa, mapq_sa, nm_sa = read.get_tag(
            'SA').split(',')
        if int(mapq_sa) < mapq_cutoff:
            return 'NA', []
    else:
        return 'NA', []
    read_mode, read_match, read_clipped = detect_softclip_mode(
        read.cigarstring)
    sa_mode, sa_match, sa_clipped = detect_softclip_mode(cigar_sa)
    if read_mode == 3 or sa_mode == 3:
        return 'NA', []
    target_start = 0
    target_end = 0
    newcigar = []
    newpos = -1
    if chrm == chr_sa:
        if strand_sa == strand:  # deletion, insertion, duplication
            if read_mode == 2 and sa_mode == 1:
                target_start = int(pos_sa) - 1
                target_end = read.aend
                if sa_clipped != 0:
                    newcigar.append([4, sa_clipped])
                newcigar.append([0, sa_match])
                target_offset = target_end - target_start
                query_offset = read.rlen - sa_clipped - read_clipped
                indel_size = query_offset - target_offset
                newpos = target_start
                if indel_size == 0:
                    return 'NA', []
                elif indel_size < 0:  # deletion
                    if abs(indel_size) > max_del_len:
                        return 'NA', []
                    newcigar.append([2, -1*indel_size])
                    newcigar.append([0, query_offset-sa_match])
                elif indel_size >= query_offset:
                    return 'NA', []
                else:  # insertion
                    newcigar.append([1, indel_size])
                    last_match_size = read.aend - (int(pos_sa)-1+sa_match)
                    if last_match_size <= 0:
                        for each in newcigar:
                            if each[0] == 0:
                                each[1] = each[1]+last_match_size
                                break
                    else:
                        newcigar.append([0, last_match_size])
                if read_clipped != 0:
                    newcigar.append([4, read_clipped])
            elif read_mode == 1 and sa_mode == 2:
                target_start = read.pos
                target_end = int(pos_sa) - 1 + sa_match
                if read_clipped != 0:
                    newcigar.append([4, read_clipped])
                newcigar.append([0, read_match])
                target_offset = target_end - target_start
                query_offset = read.rlen - sa_clipped - read_clipped
                indel_size = query_offset - target_offset
                newpos = target_start
                if indel_size == 0:
                    return 'NA', []
                elif indel_size < 0:
                    if abs(indel_size) > max_del_len:
                        return 'NA', []
                    newcigar.append([2, -1*indel_size])
                    newcigar.append([0, query_offset-read_match])
                elif indel_size >= query_offset:
                    return 'NA', []
                else:
                    newcigar.append([1, indel_size])
                    last_match_size = int(pos_sa)-1+sa_match - read.aend
                    if last_match_size <= 0:
                        for each in newcigar:
                            if each[0] == 0:
                                each[1] = each[1]+last_match_size
                                break
                    else:
                        newcigar.append([0, last_match_size])
                if sa_clipped != 0:
                    newcigar.append([4, sa_clipped])
            else:
                return 'NA', []
        else:
            return 'NA', []
    else:
        return 'NA', []
    return list(map(tuple, newcigar)), newpos


def softclipping_realignment(mapq_cutoff, max_del_len, input, output, ref_genome, gtf, splice_bin):
    bwa_bam = pysam.AlignmentFile(input, 'rb')
    output_bam = pysam.AlignmentFile('{}.temp.bam'.format(output), 'wb', template=bwa_bam)

    # for RNAseq
    splice_motif = ['GTAG', 'CTAC', 'GCAG', 'CTGC', 'ATAC', 'GTAT']
    try:
        fastafile = pysam.Fastafile(ref_genome)
    except IOError as e:
        print('read reference genome '+ref_genome+' error!', e)
        sys.exit(1)
    try:
        cvg = extract_splice_sites(gtf, splice_bin)
    except IOError as e:
        print('read GTF file '+gtf+' error!', e)
        sys.exit(1)

    try:
        for read in bwa_bam.fetch(until_eof=True):
            if read.mapq >= mapq_cutoff and not read.is_secondary and not read.has_tag('XA'):
                chrm = bwa_bam.getrname(read.rname)
                newcigar, newpos = detect_sv_from_cigar(
                    chrm, read, mapq_cutoff, max_del_len)
                if newcigar != 'NA' and newcigar != read.cigar:
                    old_cigarstring, old_cigar, old_pos = read.cigarstring, read.cigar, read.pos
                    read.cigar, read.pos = newcigar, newpos
                    if 'D' in read.cigarstring:
                        if read.is_reverse:
                            strand = '-'
                        else:
                            strand = '+'
                        junc_start, junc_end = read.blocks[0][1], read.blocks[1][0]
                        htpos1 = HTSeq.GenomicPosition(chrm, junc_start, '.')
                        htpos2 = HTSeq.GenomicPosition(chrm, junc_end, '.')
                        if cvg[htpos1] > 0 or cvg[htpos2] > 0:
                            read.cigar, read.pos = old_cigar, old_pos
                            read.setTag('JM', 'GTF')
                            output_bam.write(read)
                            continue
                        m1 = fastafile.fetch(chrm, junc_start, junc_start+2)
                        m2 = fastafile.fetch(chrm, junc_end-2, junc_end)
                        motif = m1.upper()+m2.upper()
                        if motif in splice_motif:
                            read.cigar, read.pos = old_cigar, old_pos
                            read.setTag('JM', motif)
                            output_bam.write(read)
                            continue
                    read.setTag('OA', str(old_pos+1)+','+old_cigarstring)
            output_bam.write(read)
    except ValueError as e:
        print('Bam index file is not found!', e, file=sys.stderr)
        sys.exit(1)
    bwa_bam.close()
    output_bam.close()
    try:
        subprocess.check_call(
            "samtools sort {0}.temp.bam -o {0}".format(output), shell=True, stdin=subprocess.PIPE)
    except subprocess.CalledProcessError as e:
        print('Execution failed for samtools:', e, file=sys.stderr)
        sys.exit(1)

    subprocess.check_call("samtools index {}".format(output), shell=True)

def parse_args():
    class LengthAction(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            if values <= 0:
                parser.error("Minimum length for {0} is 1".format(option_string))
            setattr(namespace, self.dest, values)
    parser = argparse.ArgumentParser(
        description = "%(prog)s -i input_bam_file -o output_bam_file [opts]",
        epilog=textwrap.dedent('''transIndel:a splice-aware algorithm that detects the mid-sized insertions
            and large deletions using chimeric alignments from DNA-seq or RNA-seq data.
            Author: Rendong Yang <yang4414@umn.edu>, MSI, University of Minnesota, 2014'''))
    parser.add_argument('-i', '--input', action='store', dest='input', help="Input BAM file", required=True)
    parser.add_argument('-o', '--output', action='store', dest='output', help="Output BAM file", required=True)
    parser.add_argument('-r', '--ref', action='store', dest='ref', help="reference genome used for analyzing RNA-seq data", required=True)
    parser.add_argument('-g', '--gtf', action='store', dest='gtf', help="gene annotatino file used for analyzing RNA-seq data", required=True)
    parser.add_argument('-s', '--splice_bin', action='store', dest='splice_bin', type=int, help="splice site half bin size (default: %(default)s)", default=20)
    parser.add_argument('-m', '--mapq', action='store', dest='mapq', type=int, help="minimal MAPQ of read from BAM file for supporting Indel (default: %(default)s)", default=15)
    parser.add_argument('-l', '--length', action=LengthAction, dest='length', type=int, help="Maximum deletion length to be detected (default: %(default)s)", default=1000000)
    parser.add_argument('-v', '--version', action='version', version='%(prog)s {}'.format(__version__))
    return parser


def external_tool_checking():
    """checking dependencies are installed"""
    software = ['samtools']
    cmd = "which"
    for each in software:
        try:
            path = subprocess.check_output(
                [cmd, each], stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
        except subprocess.CalledProcessError:
            print("Checking for '" + each + \
                "': ERROR - could not find '" + each + "'", file=sys.stderr)
            print("Exiting.", file=sys.stderr)
            sys.exit(0)
        print("Checking for '" + each + "': found " + path.decode('utf8'))
    return path.decode('utf8').rstrip()


def version_checking(samtools_path):
    '''HTSlib-based Samtools is needed!
    '''
    try:
        response = subprocess.check_output(
            '{} --version-only'.format(samtools_path), shell=True, stderr=subprocess.STDOUT, stdin=subprocess.PIPE)
    except subprocess.CalledProcessError:
        print("HTSlib-based Samtools is not available, please install latest Samtools at http://www.htslib.org/", file=sys.stderr)
        print("Exiting.", file=sys.stderr)
        sys.exit(0)


def main():
    if sys.version_info < (3,6):
        sys.exit('Sorry, this code need Python 3.6 or higher. Please update. Aborting...')
    parser = parse_args()
    if len(sys.argv[1:]) < 1:
        parser.print_help()
        sys.exit(1)
    else:
        options = parser.parse_args()

    print('transIndel_RNA build starts running: ' + time.strftime("%Y-%m-%d %H:%M:%S"))
    start = time.time()

    # check external tools used
    path = external_tool_checking()
    version_checking(path)

    # CIGAR string refinement or add SV tag
    softclipping_realignment(options.mapq, options.length,
                             options.input, options.output, options.ref, options.gtf, options.splice_bin)

    remove('{}.temp.bam'.format(options.output))
    print("transIndel_RNA build running done: " + time.strftime("%Y-%m-%d %H:%M:%S"))
    end = time.time()
    print('transIndel_RNA build takes ' + str(end - start) + ' seconds.')


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        sys.stderr.write("User interrupt me ^_^ \n")
        sys.exit(1)
