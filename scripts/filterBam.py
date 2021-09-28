#!/bin/python
# simply and quickly filter reads

import pysam
import sys, re, random, os
import argparse
import string

parser = argparse.ArgumentParser()
parser.add_argument('-b',help='Unfiltered bam file')
parser.add_argument('-o',type=str,help='Filtered but not deduped bam.')
args = parser.parse_args()


suffix = ''.join(map(str,random.sample(string.ascii_letters+string.digits,10))) + '.bam'
tmp_file = str(args.o) + '.' + suffix
abs_path = os.path.join(os.getcwd(), tmp_file)

bam = pysam.AlignmentFile(args.b,'rb')
out = pysam.AlignmentFile(tmp_file,'wb',template=bam)

for read in bam:
	if read.is_unmapped:
		continue
	if read.is_secondary:
		continue
	if read.mapping_quality < 60:
		continue
	if read.is_qcfail:
		continue
	if read.mate_is_unmapped:
		continue

	out.write(read)

bam.close()
out.close()

# sort bam
pysam.sort("-o", args.o, tmp_file)
pysam.index(args.o)

# remove tmp_file
os.remove(abs_path)
