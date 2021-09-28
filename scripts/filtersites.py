#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys
import argparse, gzip
import pandas as pd
from collections import defaultdict

parser = argparse.ArgumentParser()
parser.add_argument('-i', required=True, help='File generated from plreads.py.')
parser.add_argument('-o',required=True, help='File to save filtered results.')
parser.add_argument('-bl', required=False, help='Blacklist regions.')
args = parser.parse_args()


# make blacklists ref
bl_ref = dict()

if args.bl:
	bl = args.bl
	bl_dict = defaultdict(list)
	with gzip.open(bl) as fi:
		for line in fi:
			line = line.decode().strip().split('\t')
			chrom = line[0]
			start = int(line[1])
			end = int(line[2])
			iv = tuple([start, end])
			bl_dict[chrom].append(iv)

	bl_ref = dict()
	for k,v in bl_dict.items():
		bl_ref[k] = pd.arrays.IntervalArray.from_tuples(v)



fo = open(args.o, 'w')
with open(args.i) as fi:
	header = fi.readline().rstrip()
	print(header, file=fo)
	header = header.split('\t')

	for line in fi:
		line = line.rstrip().split('\t')
		chrom = line[header.index('chr')]
		position = int(line[header.index('position')])
		# check if position located in blacklists
		if chrom in bl_ref.keys():
			res = bl_ref[chrom].contains(position)
			if True in res:
				continue

		r2u = float(line[header.index('r2u')])
		r2d = float(line[header.index('r2d')])
		r2c = float(line[header.index('r2c')])
		Treat_total = int(float(line[header.index('Treat_total')]))

		if r2u>=2 and r2d>=2 and r2c>=2 and Treat_total >=5:
			print('\t'.join(map(str, line)), file=fo)

		else:
			continue

fo.close()
