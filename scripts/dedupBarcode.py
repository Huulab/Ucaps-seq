#!/bin/python
# remove duplicates by barcode

import pysam, sys, os
import argparse, re, string, random
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--bam',help='Bamfiles generated from filterBam.py')
parser.add_argument('--region',help='Targeted amplicon PCR genomic region. e.g. chr1:100-200')
parser.add_argument('--primer',nargs=2, help='Primers used for PCR, [ forward reverse ]')
parser.add_argument('--out', help='File name to save deduped reads with bam foramt.')
args = parser.parse_args()

# parse region
pattern = re.compile(r'(chr[0-9XYM]\d?):(\d*)-(\d*)')
pmatch = pattern.search(args.region)
chrom = pmatch.group(1)
start = int(pmatch.group(2))
end = int(pmatch.group(3))


#
fwd_p, rev_p =args.primer
bamfile = pysam.AlignmentFile(args.bam,'rb')
records = bamfile.fetch(chrom,start,end)

# prepare out file
suffix = ''.join(map(str,random.sample(string.ascii_letters+string.digits,12))) + '.bam'
tmp_file = str(args.out) + '.' + suffix
abs_path = os.path.join(os.getcwd(), tmp_file)


# validating file name
try:
	assert os.path.exists(abs_path) == False
except AssertionError:
	suffix = ''.join(map(str,random.sample(string.ascii_letters+string.digits,10))) + '.bam'
	tmp_file = str(args.o) + '.' + suffix
	abs_path = os.path.join(os.getcwd(), tmp_file)

outbam = pysam.AlignmentFile(tmp_file,'wb',template=bamfile)


total = 0
dup = 0
poses = defaultdict(int)
r1_barcodes = set()
r2_barcodes = set()

for read in records:
	total += 1
	raw_seq = read.get_forward_sequence()
	cigar = read.cigarstring
	barcode = read.query_name.rstrip().split(':')[-1]

	if read.is_reverse:
		strand = '-'
		using_primer = rev_p

	else:
		strand = '+'
		using_primer = fwd_p

	if using_primer in raw_seq:
		idx = raw_seq.index(using_primer)
		poses[idx] += 1

		if read.is_read1 and read.is_duplicate and barcode in r1_barcodes:
			dup += 1
			continue

		if read.is_read2 and read.is_duplicate and barcode in r2_barcodes:
			dup += 1
			continue
			
		else:
			if read.is_read1:
				r1_barcodes.add(barcode)
			if read.is_read2:
				r2_barcodes.add(barcode)

			outbam.write(read)
			continue
	else:
		if read.is_read1 and barcode in r1_barcodes:
			dup += 1
			continue

		if read.is_read2 and barcode in r2_barcodes:
			dup += 1
			continue

		else:
			continue



bamfile.close()
outbam.close()


# sort bam
pysam.sort("-o", args.out, tmp_file)
pysam.index(args.out)
os.remove(abs_path)


args_info = 'BAM:{}    REGION:{}    PRIMER:{}    OUTFILE:{}    MISMATCHOUT:{}'.format(str(args.bam), str(args.region), str(fwd_p+'\t'+rev_p), str(args.out), str(args.mismatchout))
print(args_info)
print('In selected region:')
print('Of total %s, duplicates:%s' % (total, dup))
print('Primer position:\n',poses)
print("")
