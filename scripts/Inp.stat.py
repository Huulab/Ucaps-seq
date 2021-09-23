#!/usr/bin/env python
# chrdU 14S99M37S
# chrdT 29S99M22S

import argparse, pysam, sys, re, json
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True, help="bam file with bai index")
parser.add_argument("--out", required=True, help="file name to save statistical results")
args = parser.parse_args()

pat = re.compile(r'^(\d{2})S99M')  # meet 99M and >10 S
total = 0 # refers to cleaned Read 1 counts
du_p,du_m,dt_p,dt_m = 0,0,0,0  # refers to strand which reads mapping to 


chroms = ['chrdU','chrdT']
bamfile = pysam.AlignmentFile(args.bam,'rb')


for read in bamfile.fetch('chrdU'):
	if read.is_read2:
		continue
	if read.mapping_quality < 60:
		continue
	if read.is_unmapped:
		continue
	if read.mate_is_unmapped:
		continue
	if read.is_secondary:
		continue
	if not read.is_proper_pair:
		continue

	cigar = str(read.cigarstring)
	m = pat.match(cigar)
	if read.is_reverse:
		if m:
			total += 1
			du_m += 1
		else:
			continue
	else:
		if m:
			total += 1
			du_p += 1
		else:
			continue


for read in bamfile.fetch('chrdT'):
	if read.is_read2:
		continue
	if read.mapping_quality < 60:
		continue
	if read.is_unmapped:
		continue
	if read.mate_is_unmapped:
		continue
	if read.is_secondary:
		continue
	if not read.is_proper_pair:
		continue

	cigar = str(read.cigarstring)
	m = pat.match(cigar)
	if read.is_reverse:
		if m:
			total += 1
			dt_m += 1
		else:
			continue
	else:
		if m:
			total += 1
			dt_p += 1
		else:
			continue

bamfile.close()

with open(args.out,'w') as fo:
	print('Of Total Reads:\t%d' % total,file=fo)
	print('dU-forward:%d\tdU-reverse:%d' % (du_p,du_m),file=fo)
	print('dT-forward:%d\tdT-reverse:%d' % (dt_p,dt_m),file=fo)
