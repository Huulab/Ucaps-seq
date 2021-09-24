#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse,sys,gzip
import tabix, py2bit

parser = argparse.ArgumentParser()

parser.add_argument('-i',required=True,type=str,nargs=2,help='Treat.bed Input.bed')
parser.add_argument('-b',required=True,help='RT_HeLaS3_hg19.10kb.bedgraph')
parser.add_argument('-g',required=True,help='hg19.2bit')
parser.add_argument('-o',help='File name to save results')

args = parser.parse_args()


def get_base_content(chrom,start,end,tbit,fraction=True,base='G'):
	valid_bases = ['A','C','G','T']
	base = base.upper()
	if base not in valid_bases:
		raise ValueError("Invalid base type. Please input A, C, G or T [ignore case].")
		exit()
	bases = tbit.bases(chrom,start,end,fraction=False)
	if end > tbit.chroms(chrom):
		end = tbit.chroms(chrom)
	if sum(bases.values()) < 0.95 * (start - end):
		raise Exception("WARNING: too many NNNs present in {}:{}-{}".format(chrom, start, end))
		return None
	if fraction:
		return (bases[base]) / float(end - start)
	return bases[base]


def get_cov(records,start,end):

	p,m = 0,0
	for read in records:
		read_strand = read[-1]
		read_hit = int(read[1])
		if read_hit<start or read_hit>end:
			continue

		if read_strand == '+':
			p += 1
		if read_strand == '-':
			m += 1
		else:
			continue
			
	return p,m


def main():

	treat,inp = args.i
	tbed = tabix.open(treat)
	ibed = tabix.open(inp)
	tbit = py2bit.open(args.g)

	out = open(args.o,'w')

	with open(args.b) as fi:
		for line in fi:
			line = line.rstrip().split('\t')
			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			score = float(line[3])
			trecords = tbed.query(chrom,start,end)
			irecords = ibed.query(chrom,start,end)
			tplus,tminus = get_cov(trecords,start,end)
			iplus,iminus = get_cov(irecords,start,end)
			base_T = get_base_content(chrom,start,end,tbit,fraction=False,base='T')
			base_A = get_base_content(chrom,start,end,tbit,fraction=False,base='A')

			print('\t'.join(map(str,['\t'.join(map(str,line)),tplus,tminus,iplus,iminus,base_T,base_A])), file=out)

	tbit.close()
	out.close()


if __name__ == '__main__':
	main()
