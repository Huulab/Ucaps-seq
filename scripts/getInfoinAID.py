#!/usr/bin/env python
# -*- coding:utf-8 -*-

import pysam,argparse
import re
from pyfaidx import Fasta

parser = argparse.ArgumentParser()
parser.add_argument('--ibam',help='Bam file with bai index')
parser.add_argument('--region', type=str, nargs='+', help='Genome region in `chr12:100-200` format.')
parser.add_argument('--genome',help='Reference genome.[mm10 fasta]')
parser.add_argument('--out',help='File name to save results.')
args = parser.parse_args()


def fetch_sequence(fasta, chrom, start, end):
    try:
        sequence = fasta[chrom][start:end]
    except KeyError:
        sys.stderr.write("warning: {name} not found in file\n".format(**locals()))
        return
    except ValueError as ve:
        sys.stderr.write(str(ve))
        return
    return sequence


def check_C(seq):
    if seq[2] == 'C':
        return True
    else:
        return False


motif_WRC = ['AAC','AGC','TAC','TGC']

fasta = Fasta(args.genome)
bamfile = pysam.AlignmentFile(args.ibam,'rb')
out = open(args.out,'w')

total_mapped = bamfile.mapped
pattern = re.compile(r'^chr12:(\d+)-(\d+)$')



for line in args.region:

    base_counts = {'A':0,'C':0,'G':0,'T':0}

    line = line.rstrip()
    match = pattern.match(line)
    start = int(match.group(1))
    end = int(match.group(2))

    raw_count = 0
    filter_count = 0
    C_count = 0
    agct_count = 0
    wrc_count = 0

    for read in bamfile.fetch('chr12', start, end):

        raw_count += 1

        if read.is_duplicate:
            continue
        if read.is_read2:
            continue
        if read.is_unmapped:
            continue
        if read.mate_is_unmapped:
            continue
        if read.mapping_quality < 60:
            continue
        if read.is_secondary:
            continue
        if not read.is_proper_pair:
            continue

        if read.is_reverse:
            if re.match(".*\d+S$", read.cigarstring):
                continue

            r_start = read.reference_end
            if r_start <start or r_start > end:
                continue

            sequence = fetch_sequence(fasta,'chr12',r_start-1,r_start+3)
            if sequence == None:
                continue

            sequence = sequence.reverse.complement.seq.upper()
            if sequence.count('N'):
                continue

            filter_count += 1

        else:
            if re.match("^\d+S.*", read.cigarstring):
                continue

            r_start = read.reference_start
            if r_start <start or r_start > end:
                continue

            sequence = fetch_sequence(fasta,'chr12',r_start-3,r_start+1)
            if sequence == None:
                continue

            sequence = sequence.seq.upper()
            if sequence.count('N'):
                continue

            filter_count += 1

        if check_C(sequence):
            C_count += 1

        base_counts[sequence[2]] += 1

        if sequence[:3] in motif_WRC and sequence != 'AGCT':
            wrc_count += 1

        if sequence == 'AGCT':
            agct_count += 1

        inmotif = wrc_count + agct_count

        for k,v in base_counts.items():
            print('{}\t{:.0f}'.format(k,v),file=out)

        print('%s\t%s\t%s\t%s\t%s\t%s\t%s' % (total_mapped,raw_count,filter_count,C_count,inmotif,wrc_count,agct_count), file=out)


out.close()
bamfile.close()
