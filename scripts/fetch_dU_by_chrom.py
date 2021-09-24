#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse, pysam, sys, re
from pyfaidx import Fasta, FetchError


def get_arguments():

	parser = argparse.ArgumentParser()
	parser.add_argument("--ref_fa", required=True, help="fasta file with fai index")
	parser.add_argument("--bam", required=True, help="input bam file with bai index")
	parser.add_argument("--out_file", required=True, help="output file with bed format")
	parser.add_argument("--chrom", required=True, help="chrom[chr1,...]")
	
	return parser


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


def main():

	args = get_arguments().parse_args()
	ref_fa = args.ref_fa
	bam = args.bam
	out_file = args.out_file
	chrom = args.chrom

	fasta = Fasta(ref_fa)
	samfile = pysam.AlignmentFile(bam, "rb")
	out = open(out_file, "w")

	# fetch and filter reads
	for read in samfile.fetch(chrom):	
		if read.is_read2:
			continue
		if read.is_unmapped:
			continue
		if read.mate_is_unmapped:
			continue
		if not read.is_proper_pair:
			continue
		if read.is_secondary:
			continue
		if read.mapping_quality < 60:
			continue

		# for negtive strand
		# dU located in read end + 1
		if read.is_reverse:
			# removing reads contains customed adapter
			if re.match(".*\d+S$", read.cigarstring):
				continue

			start = read.reference_end
			end = read.reference_end + 1
			strand = "-"

			sequence = fetch_sequence(fasta, chrom, start, end)
			if sequence == None:
				continue
			seq_up = sequence.reverse.complement.seq.upper()
			score = 1000 - mapq

		# for positive strand
		# dU located in read start - 1
		else:
			if re.match("^\d+S.*", read.cigarstring):
				continue

			start = read.reference_start - 1
			end = read.reference_start
			strand = "+"

			sequence = fetch_sequence(fasta, chrom, start, end)
			if sequence == None:
				continue
			seq_up = sequence.seq.upper()
			score = mapq

		if seq_up == "N":
			continue

		print("\t".join(map(str, [chrom, start, end, seq_up, score, strand])),file=out)


	out.close()
	samfile.close()
	return

if __name__ == '__main__':
	main()

