#!/usr/bin/env python
# chrdU 14S99M37S
# chrdT 29S99M22S

import argparse, pysam, sys, re, json
from collections import defaultdict


parser = argparse.ArgumentParser()
parser.add_argument("--bam", required=True, help="bam file with bai index")
parser.add_argument("--out", required=True, help="file name to save statistical results")
parser.add_argument('--json',help='statistic cigars for those reads in keeps.')
args = parser.parse_args()


pat = re.compile(r'^(\d{2})S99M')  # meet 99M and >10 S
total = 0 # refers to cleaned Read 1 counts
specific, non_specific = 0, 0 # refers to read counts mapping to dU or dT
du_p,du_m,dt_p,dt_m = 0,0,0,0  # refers to strand which reads mapping to 
to_keep = list(range(11,18))  # limits dU range

fp = 0 # false positive: refers to reads counted in dU but not real dU
exceed = 0 # mapped to dU forward and cigarstring startswith 15/16/17S


def check_T(seq,sv,fp,exceed):
	if sv <= 14:
		return fp,exceed
	else:
		exceed += 1
		idx = sv - 15
		base = seq[idx]
		if base == 'T':
			fp += 1
		return fp,exceed


fwd_collect = {}
for i in to_keep:
	fwd_collect[str(i)+'S'] = defaultdict(int)


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
	if read.is_reverse:
		m = pat.match(cigar)
		if m:
			total += 1
			du_m += 1
			non_specific += 1
		else:
			continue
	else:
		sequence = read.get_forward_sequence()
		m = pat.match(cigar)
		if m:
			total += 1
			du_p += 1
			soft_clip = int(m.group(1))
			if soft_clip in to_keep:
				specific += 1
				fwd_collect[str(soft_clip)+'S'][cigar] += 1
				fp,exceed = check_T(sequence,soft_clip,fp,exceed)
			else:
				non_specific += 1
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
			non_specific += 1
		else:
			continue
	else:
		if m:
			total += 1
			dt_p += 1
			non_specific += 1
		else:
			continue


bamfile.close()


with open(args.out,'w') as fo:
	print('Of Total Reads:\t%d' % total,file=fo)
	print('Specific:%d\tnon-Specific:%d' % (specific,non_specific),file=fo)
	print('dU-forward:%d\tdU-reverse:%d' % (du_p,du_m),file=fo)
	print('dT-forward:%d\tdT-reverse:%d' % (dt_p,dt_m),file=fo)
	print('Within Specific:\n\tExceed:%d\tTrue positive:%d' % (exceed,fp),file=fo)

with open(args.json,'w') as f:
	json.dump(fwd_collect,f,indent=4)
