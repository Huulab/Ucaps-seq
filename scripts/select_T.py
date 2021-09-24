#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys,gzip

bed = sys.argv[1]

if str(bed).endswith('.gz'):
	with gzip.open(bed) as f:
		for line in f:
			line = line.decode().rstrip().split('\t')
			sequence = str(line[3])
			if not sequence == 'T':
				continue
			else:
				print('\t'.join(map(str,line)))


else:
	with open(bed) as f:
		for line in f:
			line = line.rstrip().split('\t')
			sequence = str(line[3])
			if not sequence == 'T':
				continue
			else:
				print('\t'.join(map(str,line)))
