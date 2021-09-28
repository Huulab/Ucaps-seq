#!/bin/python
# merge mutant infomation

import pysam, sys, argparse
import os
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i',help='Mutant info file generated from getMutInfo.py')
parser.add_argument('--bam', help='Bam file used in getMutInfo.py')
parser.add_argument('-o',help='Merged mutant info')
args = parser.parse_args()


def get_chroms(infile):

	cmd = "cat %s |awk '{print $2}' |sort | uniq" % infile
	chroms = list(os.popen(cmd).readlines())
	chroms = [x.strip() for x in chroms]

	return chroms


def get_coverage(bamfile, chrom, pos, strand):

	total_count = bamfile.count(chrom, pos, pos+1, read_callback='nofilter')
	records = bamfile.fetch(chrom,pos,pos+1)
	plus, minus = 0, 0

	for read in records:
		if read.is_reverse:
			minus += 1
		else:
			plus += 1

	try:
		assert total_count == plus+minus
	except AssertionError:
		print('\t'.join(map(str, [chrom, pos])))
		print('total: %d' % total_count)
		print('plus:%d\tminus:%d' % (plus, minus))
		sys.exit()

	if strand == '+':
		return plus, total_count
	else:
		return minus, total_count


def main():

	mut_file = args.i
	chroms = get_chroms(mut_file)

	res = {}
	for strand in ['+','-']:
		res[strand] = {}
		for chrom in chroms:
			res[strand][chrom] = defaultdict(int)


	with open(mut_file) as fi:
		for line in fi:
			line = line.rstrip().split('\t')
			sub_chrom = line[1]
			sub_pos = line[2] # 0-based
			ref_base = line[3]
			mut_base = line[4]
			sub_strand = line[5]

			mut_info = '_'.join(map(str, [sub_pos, ref_base, mut_base]))
			res[sub_strand][sub_chrom][mut_info] += 1


	bamfile = pysam.AlignmentFile(args.bam)
	with open(str(args.o)+'.tmp', 'w') as fo:
		print('\t'.join(map(str, ['chrom','0_based_pos','ref','mut','count','strand','strand_counts','all_counts'])), file = fo)

		for strand in ['+','-']:
			for chrom in chroms:
				tmp = res[strand][chrom]
				mut_list = tmp.keys()
				for mut_type in mut_list:
					count = tmp.get(mut_type)
					dezip = list(mut_type.split('_'))
					strand_total, total_count = get_coverage(bamfile, chrom, int(dezip[0]), strand)
					dezip.insert(0,chrom)
					dezip.append(count)
					dezip.append(strand)
					dezip.append(strand_total)
					dezip.append(total_count)
					print('\t'.join(map(str, dezip)), file = fo)

	bamfile.close()

	cmd_1 = 'sort -k5,5nr %s > %s' % (str(args.o)+'.tmp', str(args.o))
	cmd_2 = 'rm %s' % str(args.o)+'.tmp'
	os.system(cmd_1)
	os.system(cmd_2)


if __name__ == '__main__':
	main()
