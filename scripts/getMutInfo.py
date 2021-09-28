#!/bin/python

import pysam, sys, argparse
import os,pysnooper
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('--bam',help='Bam file generated from dedupBarcode.py')
parser.add_argument('--out',help='Mutant info')
parser.add_argument('--target', required=False, nargs='+', help='keep specific mutant type. e.g. CT')
args = parser.parse_args()


def transfer_base(a,b):
	base_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A', 'N':'N'}
	new_a = base_dict.get(a)
	new_b = base_dict.get(b)

	return new_a, new_b

def process_reads(read, keep_mutant):

	mut_info = []

	# strand
	strand = '+'
	if read.is_reverse:
		strand = '-'

	# chromosome
	sub_chrom = read.reference_name

	pairs = read.get_aligned_pairs(with_seq=True)
	seq = read.query_sequence
	q_name = read.query_name

	#
	for match_pair in pairs:
		if match_pair[0] == None or match_pair[1] == None or match_pair[2] == None:
			continue

		if match_pair[2].islower():
			q_idx = int(match_pair[0])
			ref_pos = int(match_pair[1]) # 0-based coordinates
			ref_base = str(match_pair[2]).upper()
			query_base = seq[q_idx].upper()

			if strand == '-':
				ref_base, query_base = transfer_base(ref_base, query_base)

			sub_type = str(ref_base) + str(query_base)
			if len(keep_mutant)!=0 and sub_type not in keep_mutant:
				continue

			sub_mut_info = '_'.join(map(str, [q_name, sub_chrom, ref_pos, ref_base, query_base, strand]))
			mut_info.append(sub_mut_info)

		else:
			continue

	return mut_info



def main():

	keep_mutant = []
	if args.target:
		for mut in args.target:
			if len(mut) != 2:
				sys.stderr.write('%s\t is not a valid format.' % mut)
			else:
				keep_mutant.append(mut)


	# create dicts for data
	bamfile = pysam.AlignmentFile(args.bam)
	merge_info = []

	while True:
		try:
			read = next(bamfile)
			mut_info = process_reads(read, keep_mutant)
			merge_info.extend(mut_info)

		except StopIteration:
			break

	bamfile.close()

	with open(args.out,'w') as fo:
		for item in merge_info:
			out_line = item.split('_')
			print('\t'.join(map(str, out_line)), file=fo)


if __name__ == '__main__':
	main()			
