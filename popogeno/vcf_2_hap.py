#! /usr/bin/python

import sys
import argparse
import os.path
import re
import vcf
from Bio.Alphabet import generic_dna
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO

#define global variables
formats = []
position = []
inddict = {}
pops_list = []
traitstring_pops = {}
seqdict={}
seqrecords = []
w_count = 1

###DEFINING FUNCTIONS
def add_traits_block(p_list, i_dict, nexhandle):
	"add a traits block to the end of a (nexus) file"
	traitstring_pops = {}
	for i in range(len(p_list)):
		traitstring = []
		for j in range(len(p_list)):
			if j == i:
				traitstring.append("1")
			else:
				traitstring.append("0")
		traitstring_pops[p_list[i]] = ",".join(traitstring)

	nexhandle.write("\nBEGIN TRAITS;\n")
	nexhandle.write(" Dimensions NTRAITS=%i;\n" %len(p_list))
	nexhandle.write(" Format labels=yes missing=? separator=Comma;\n")
	nexhandle.write(" Traitlabels %s;\n" %" ".join(p_list))
	nexhandle.write(" Matrix\n")
	for ind in i_dict.keys():
		nexhandle.write("%s_0 %s\n" %(ind, traitstring_pops[i_dict[ind]]))
		nexhandle.write("%s_1 %s\n" %(ind, traitstring_pops[i_dict[ind]]))
	nexhandle.write(";\n\nEND;")


###
parser = argparse.ArgumentParser(description='produce fasta/nexus file formatted haplotypes from vcf input file')
parser.add_argument('in_vcf', metavar='<VCF file>', help='vcf file to be translated into haplotypes (gzipped possible)')
parser.add_argument('-p', '--popmap', help='Tab delimited text file to assign individuals to populations. col1 = population ID, col2 = individual (as in vcf).', metavar="<FILE>", type=str, action="store")
parser.add_argument('-n', '--nexus', action='store_true', help='output haplotypes in nexus format')
parser.add_argument('-f', '--fasta', action='store_true', help='output haplotypes in fasta format')
parser.add_argument('-w', '--window', action='store', metavar='<int>', help='window size in kb', type=int)
parser.add_argument('-s', '--step', action='store', metavar='<int>', help='stepsize in kb', type=int)
parser.add_argument('--minlength', action='store', metavar='<int>', help='minimum length of haplotype to be written', type=int, default=1)
parser.add_argument('-o', '--out_prefix', action='store', metavar='<string>', help='prefix for output files')
args = parser.parse_args()

#check input
if len(sys.argv) < 2:   #if the script is called without any arguments display the usage
    parser.print_usage()
    sys.exit(1)

if not args.fasta and not args.nexus:
	print "Please specify desired output format\n"
	sys.exit()
else:
	if args.fasta:
		formats.append('fasta')
	if args.nexus:
		formats.append('nexus')

if not os.path.isfile(args.in_vcf):
	print "Please provide path to valid vcf file"
	sys.exit(0)
else:
	vcf = vcf.Reader(open(args.in_vcf, 'r'))

if not args.out_prefix:
	args.out_prefix='out'

if (args.window and not args.step) or (not args.window and args.step):
	print "both window and stepsize needs to be specified\n"
	sys.exit(0)

#parse populationmap
if args.popmap:
        if not os.path.isfile(args.popmap):
                print "The populationmap provided at: %s is not a valid file" %args.popmap
                sys.exit(0)
        else:
                pops = open(args.popmap ,'r')
                for pop in pops.readlines():
			pop = pop.strip()
			ind = pop.split("\t")[0]
			pop = pop.split("\t")[1]
			if not ind and pop:
				print "The populationmap is expected to have at least 2 tab delimited columns"
				sys.exit(0)

			inddict[ind] = pop
			pops_list.append(pop)
	
			if not seqdict.has_key(ind):
				seqdict[ind]=["",""]
		
		pops_list = sorted(set(pops_list))
else:
	print "please provide a populationmap"

#parse vcf file
for record in vcf:
	allele_dict = {'0': record.REF, '1': str(record.ALT[0])}
#	print allele_dict
	position.append(record.POS)
	for sample in inddict.keys():
#		print sample
#		print record.genotype(sample)['GT']
#		print allele_dict[record.genotype(sample)['GT'].split("|")[0]]+"|"+allele_dict[record.genotype(sample)['GT'].split("|")[1]]
		seqdict[sample][0]+=allele_dict[record.genotype(sample)['GT'].split("|")[0]]
		seqdict[sample][1]+=allele_dict[record.genotype(sample)['GT'].split("|")[1]]
#		print seqdict[sample]
	
#print seqdict

#generate alignment
for ind in inddict.keys():
	for i in range(2):
#		print ">"+ind+"_"+str(i)
#		print seqdict[ind][i]
		seqrec=SeqRecord(Seq(seqdict[ind][i], generic_dna), id=ind+"_"+str(i), description=ind+"_"+str(i))
#		print seqrec
		seqrecords.append(seqrec)

align = MultipleSeqAlignment(seqrecords)

#print alignemnt in desired formats
for f in formats:
	print "writing haplotypes to "+f+"\n"
	OUT = open(args.out_prefix+'.'+f[:3],'w')
	OUT.write(align.format(f))
	if f is 'nexus':
		add_traits_block(p_list=pops_list, i_dict=inddict, nexhandle=OUT)
	OUT.close()

#produce sliding windows
if args.window and args.step:
	start = min(position)
	end = max(position)

	for w in range(start, end, args.step):
#		print "%s to %s" %(w, w+args.window)
		slice_index = []
		for i in range(len(position)):
			if position[i] >= w and position[i] <= w+args.window:
#				print "\t"+str(position[i])
				slice_index.append(i)
			elif position[i] > w+args.window:
				break

		if len(slice_index) >= args.minlength:
#			print str(min(slice_index))+":"+str(max(slice_index))
			for f in formats:
				OUT=open(args.out_prefix+'_w_'+str('%05d' %w_count)+'_'+str(w)+'-'+str(w+args.window)+'.'+f[:3], 'w')
				OUT.write(align[:, min(slice_index):max(slice_index)+1].format(f))
				if f is 'nexus':
					add_traits_block(p_list=pops_list, i_dict=inddict, nexhandle=OUT)
				OUT.close()
				w_count += 1
