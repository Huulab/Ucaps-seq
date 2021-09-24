#!/usr/bin/env python
# -*- coding:utf-8 -*-
# Usage: count_base.py xxx.bed(.gz) xxx.basecount.txt genome.fa
import sys, gzip
import pandas as pd
from pyfaidx import Fasta

def count_single_base(seq_list, index):
	base_dict = {}
	for seq in seq_list:
		base = seq[index]
		k = base_dict.get(base, 0)
		base_dict[base] = k + 1
	return base_dict

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

ref = Fasta(sys.argv[3])

# open bedfile
bed = sys.argv[1]
if not bed.endswith('gz'):
	positions = [int(line.strip().split("\t")[1]) for line in open(bed)]
	chroms = [line.rstrip().split('\t')[0] for line in open(bed)]
	strands = [line.rstrip().split('\t')[-1] for line in open(bed)]
else:
	positions = [int(line.decode().strip().split("\t")[1]) for line in gzip.open(bed)]
	chroms = [line.decode().rstrip().split('\t')[0] for line in gzip.open(bed)]
	strands = [line.decode().rstrip().split('\t')[-1] for line in gzip.open(bed)]

positions = [[x-3,x+4] for x in positions]
seq_list = []


# fetch sequence
for i in range(len(positions)):
	sub_chrom = chroms[i]
	a, b = positions[i]
	sub_strand = strands[i]
	sub_seq = fetch_sequence(ref, sub_chrom, a, b)
	seq_res = sub_seq.seq.upper()
	if sub_strand == '-':
		seq_res = sub_seq.reverse.complement.seq.upper()

	if seq_res.count('N'):
		continue

	seq_list.append(seq_res)


ref.close()
outfile = sys.argv[2]
dat = pd.DataFrame(0,index=range(7),columns=['A', 'C', 'G', 'T'])

for i in range(7):
	basestmp = count_single_base(seq_list,i)
	for k in basestmp.keys():
		dat.loc[i,k] = basestmp.get(k)

dat.to_csv(outfile,sep='\t',float_format=None,header =True,index=True, quoting=None)
