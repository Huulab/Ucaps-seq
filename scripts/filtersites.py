#!/usr/bin/env python
# -*- coding:utf-8 -*-

import sys

fo = open(sys.argv[2],'w')

with open(sys.argv[1]) as fi:
	header = fi.readline().rstrip()
	print(header, file=fo)
	header = header.split('\t')

	for line in fi:
		line = line.rstrip().split('\t')
		r2u = float(line[header.index('r2u')])
		r2d = float(line[header.index('r2d')])
		r2c = float(line[header.index('r2c')])
		Treat_total = int(float(line[header.index('Treat_total')]))

		if r2u>=2 and r2d>=2 and r2c>=2 and Treat_total >=5:
			print('\t'.join(map(str, line)), file=fo)

		else:
			continue

fo.close()
