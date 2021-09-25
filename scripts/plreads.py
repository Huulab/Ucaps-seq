#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys, argparse
import tabix, gzip
from collections import defaultdict
from pyfaidx import Fasta


def load_bed_as_dict(infile):

	keep_chroms = ['chr'+str(i) for i in range(1,23)] + ['chrX']
	res = {}
	for strand in ['+','-']:
		res[strand] = {}
		for chrom in keep_chroms:
			res[strand][chrom] = defaultdict(int)

	total = 0

	with gzip.open(infile) as fi:
		for line in fi:
			line = line.decode().strip().split('\t')
			sub_chrom = line[0]
			if not sub_chrom in keep_chroms:
				continue

			total += 1
			snu_site = int(line[2]) # 1-based coordinates
			sub_strand = line[-1]
			snu_base = line[3]

			res[sub_strand][sub_chrom][snu_site] += 1

	return total, res

def template_dict():
	keep_chroms = ['chr'+str(i) for i in range(1,23)] + ['chrX']
	res = {}
	for strand in ['+','-']:
		res[strand] = {}
		for chrom in keep_chroms:
			res[strand][chrom] = {}

	return res


def update_dict(raw, added):

	res_dict = raw.copy()

	for strand in ['+','-']:
		chroms = added[strand].keys()
		for chrom in chroms:
			tmp_dict = added[strand][chrom]
			for k,v in tmp_dict.items():
				if res_dict[strand][chrom].get(k):
					res_dict[strand][chrom][k] += v
				else:
					res_dict[strand][chrom][k] = v

	return res_dict


def cal_ratio(a,b,c,d):
	cov = a
	try:
		r2u = (a+1) / (b+1)
	except TypeError:
		r2u = a + 1

	try:
		r2d = (a+1) / (c+1)
	except TypeError:
		r2d = a + 1

	try:
		r2c = (a+1) / (d+1)
	except TypeError:
		r2c = a + 1

	return cov, r2u, r2d, r2c


def transfer_base(base):
	base_dict = {'A':'T','T':'A','G':'C','C':'G','N':'N'}
	t_base = base_dict.get(base)
	return t_base


def get_single_count(strand, chrom, pos, dict_a, dict_b):
	# dict A refers to treat samples and dict B refers to input samples
	t_count, c_count = [], []

	for t_sm in dict_a.keys():
		pos_count = dict_a[t_sm][strand][chrom].get(pos)
		if pos_count == None:
			pos_count = 0
		t_count.append(pos_count)

	for c_sm in dict_b.keys():
		pos_count = dict_b[c_sm][strand][chrom].get(pos)
		if pos_count == None:
			pos_count = 0
		c_count.append(pos_count)

	t_count.extend(c_count)
	res = t_count[:]
	return res


def main():
	parser = argparse.ArgumentParser()
	parser.add_argument('-t',nargs='+',help='Treat group bedfiles.')
	parser.add_argument('-c',nargs='+',help='INput group bedfiles.')
	parser.add_argument('-g',help='hg19.fasta')
	parser.add_argument('-o',help='output file.')
	args = parser.parse_args()

	# load and merge data
	# treat samples
	treat_title = []
	treat_samples = args.t
	treat_samplecounts = len(treat_samples)
	for i in range(1,treat_samplecounts+1):
		treat_title.append('Treat_' + str(i))

	treat_total_counts = 0
	treat_dict = template_dict()
	treat_sole_dict = {}

	for i in range(treat_samplecounts):
		sample = treat_samples[i]
		sub_count, sub_treat_dict = load_bed_as_dict(sample)
		treat_sole_dict[treat_title[i]] = sub_treat_dict.copy()
		treat_total_counts += sub_count
		treat_dict = update_dict(treat_dict, sub_treat_dict)


	# input samples
	input_title = []
	input_samples = args.c
	input_samplecounts = len(input_samples)
	for i in range(1,input_samplecounts+1):
		input_title.append('Input_' + str(i))

	input_total_counts = 0
	input_dict = template_dict()
	input_sole_dict ={}

	for i in range(input_samplecounts):
		sample = input_samples[i]
		sub_count, sub_input_dict = load_bed_as_dict(sample)
		input_sole_dict[input_title[i]] = sub_input_dict.copy()
		input_total_counts += sub_count
		input_dict = update_dict(input_dict, sub_input_dict)


	total_ratio = treat_total_counts / input_total_counts

	# load fasta file
	genome = Fasta(args.g)
	keep_chroms = ['chr'+str(i) for i in range(1,23)] + ['chrX']
	header = ['chr','position','strand', 'r2u','r2d','r2c','Treat_total','Input_total'] + treat_title + input_title
	# chrom position strand  r2u r2d r2c treat_total input_total treat_1 treat_2 ... input_1 input_2
	# 3 ratios were calculcated by the raw coverage plus 1 
	with open(args.o,'w') as fo:
		print('\t'.join(map(str, header)), file=fo)
		for strand in ['+','-']:
			for chrom in keep_chroms:
				out_dict = treat_dict[strand][chrom]
				for pos in sorted(out_dict.keys()):
					cur_base = genome[chrom][pos-1].seq.upper()
					if strand == '-':
						cur_base = transfer_base(cur_base)
					if cur_base != 'C':
						continue

					cur_cov = out_dict.get(pos)
					up_cov = out_dict.get(pos-1)
					down_cov = out_dict.get(pos+1)
					input_cov = input_dict[strand][chrom].get(pos)
					cov, r2u, r2d, r2c = cal_ratio(cur_cov, up_cov, down_cov, input_cov)
					if input_cov:
						r2c /= total_ratio # total_ratio is not considering the condition of cov+1, but it doesn't matter

					if input_cov == None:
						input_cov = 0

					single_count = get_single_count(strand, chrom, pos, treat_sole_dict, input_sole_dict)
					assert len(single_count) == treat_samplecounts+input_samplecounts
					out_line = [chrom, pos, strand, r2u, r2d, r2c, cov, input_cov] + single_count
					print('\t'.join(map(str, out_line)), file=fo)


if __name__ == '__main__':
	main()
