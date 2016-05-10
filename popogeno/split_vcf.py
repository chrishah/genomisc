#! /usr/bin/python

import sys
import argparse
import os.path
import re
import gzip



#define global variables
comments = []
calls = []
filecount=1
###DEFINING FUNCTIONS


###
parser = argparse.ArgumentParser(description='split vcf file into n sub files')
parser.add_argument('in_vcf', metavar='<VCF file>', help='input vcf file (gzipped possible)')
parser.add_argument('-n', '--number', action='store', metavar='<int>', type=int, default=100000, help='number of calls per file')
parser.add_argument('-o', '--out_prefix', action='store', metavar='<string>', default='out', help='prefix for output files, the file number will be appended to this')
args = parser.parse_args()

#check input
if len(sys.argv) < 2:   #if the script is called without any arguments display the usage
    parser.print_usage()
    sys.exit(1)

if not os.path.isfile(args.in_vcf):
	print "Please provide path to valid vcf file"
	sys.exit(0)

print "Splitting %s into chunks of %s calls -> fileprefix: %s" %(args.in_vcf, args.number, args.out_prefix) 

if args.in_vcf.endswith("gz"):
	FH = gzip.open(args.in_vcf,'rb')
else:
	FH = open(args.in_vcf,'r')

for line in FH:
	if line.startswith("#"):
#		print line
		comments.append(line)
	else:
		calls.append(line)

	if len(calls) == args.number:
		print "writing: %s-%08d.vcf" %(args.out_prefix, filecount)
		OUT = open(args.out_prefix+"-%08d.vcf" %(filecount), 'w')
		for l in comments:
			OUT.write(l)
		for c in calls:
			OUT.write(c)
		filecount+=1
		calls = []
		OUT.close()
#		OUT = open(args.out_prefix+"")

FH.close()


