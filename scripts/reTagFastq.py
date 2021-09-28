#!/usr/bin/env python
# -*- coding:utf-8 -*-

import argparse, sys, os
import gzip

parser = argparse.ArgumentParser()
parser.add_argument('-i',nargs=2,help='xxx_R1.fastq.gz xxx_R2.fastq.gz')
parser.add_argument('-p',type=str,help='prefix of tagged xxx_R1/R2.fastq.gz')
args = parser.parse_args()

'''
def check_files(a,b):
	# load files & check
	r1 = gzip.open(a)
	r2 = gzip.open(b)
	line1, line2 = 0, 0

	flag = True
	while flag:
		line = r1.readline().strip()
		if line == "" or line == None:
			flag = False
		else:
			line1 += 1

	flag = True
	while flag:
		line = r2.readline().strip()
		if line == "" or line == None:
			flag = False
		else:
			line2 += 1

	r1.close()
	r2.close()

	try:
		assert line1 == line2
	except AssertionError:
		sys.exit('2 input fastq files not same in length.')

	return line1
'''

def check_files(a,b):
	cmd_a = 'zcat %s | wc -l ' % a
	a_lines = int(os.popen(cmd_a).readline().strip())

	cmd_b = 'zcat %s | wc -l ' % b
	b_lines = int(os.popen(cmd_b).readline().strip())

	try:
		assert a_lines == b_lines
	except AssertionError:
		sys.exit('Input fastq files not same in length.')

	return b_lines


def read_by_cycle(infile):
	header = infile.readline().decode().strip()
	seq = infile.readline().decode().strip()
	anno = infile.readline().decode().strip()
	seq_qul = infile.readline().decode().strip()

	res = [header, seq, anno, seq_qul]
	return res


def main():
	file_a, file_b = args.i
	total_lines = check_files(file_a, file_b)
	abs_path = os.path.abspath(file_a)
	dir_name = os.path.dirname(abs_path)

	r1_out = gzip.open(os.path.join(dir_name,args.p+'_tagged_R1.fastq.gz'), 'wb')
	r2_out = gzip.open(os.path.join(dir_name,args.p+'_tagged_R2.fastq.gz'), 'wb')

	counter = 0
	r1 = gzip.open(file_a)
	r2 = gzip.open(file_b)

	while counter < total_lines:
		# read 4 lines per cycle
		r1_group = read_by_cycle(r1)
		r2_group = read_by_cycle(r2)
		counter += 4
		# fetch the barcode
		r1_barcode = r1_group[1][:6]
		r2_barcode = r2_group[1][:6]
		m_barcode = r1_barcode + r2_barcode
		# add barcode to header
		r1_header = r1_group[0].split(' ')
		r1_reheader = r1_header[0] + ':' + m_barcode + ' ' + r1_header[1]
		r1_group[0] = r1_reheader

		r2_header = r2_group[0].split(' ')
		r2_reheader = r2_header[0] + ':' + m_barcode + ' ' + r2_header[1]
		r2_group[0] = r2_reheader

		# write out
		r1_compress = [(x + '\n').encode() for x in r1_group]
		r2_compress = [(x + '\n').encode() for x in r2_group]
		#if counter == total_lines:
		#	r1_compress = [(x + '\n').encode() for x in r1_group[:3]] + [r1_group[3].encode()]
		#	r2_compress = [(x + '\n').encode() for x in r2_group[:3]] + [r2_group[3].encode()]

		r1_out.writelines(r1_compress)
		r2_out.writelines(r2_compress)

	r1.close()
	r2.close()
	r1_out.close()
	r2_out.close()


if __name__ == '__main__':
	main()
