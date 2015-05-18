
import sys ## Exportable and executable
import os.path
from collections import defaultdict
import random
import argparse
from argparse import RawTextHelpFormatter
import vcf

parser = argparse.ArgumentParser(description='converts vcf format to SNPfile as required by Bayenv2')

parser.add_argument('VCF', help='vcf 4.0 (gzipped is supported)')
optional_group = parser.add_argument_group('Output options', 'The parameters in this group affect the format and name of the output')
parser.add_argument('-p', '--min_prop', help='minimum proportion of individuals with data per population (default: 1.00)', metavar="<FLOAT>", type=float, action="store", default=1.0)
parser.add_argument('-n', '--min_number', help='minimum number of individuals with data per population (default: 0, which means all individuals per population)', metavar="<INT>", type=int, action="store")
parser.add_argument('-r', '--random_n', help='number of randomly sampled loci to be outputted (default = 0, which means all loci)', metavar="INT", type=int, action="store", default=0)
parser.add_argument('-m', '--popmap', help='Tab delimited text file to assign individuals to populations. col1 = population ID, col2 = individual (as in vcf).', metavar="<FILE>", type=str, action="store")

output_group = parser.add_argument_group('Output options', 'The parameters in this group affect the format and name of the output')
output_group.add_argument('-b', '--bayenv', help='output bayenv format (default)', action="store_true")
output_group.add_argument('-t', '--treemix', help='output treemix format', action="store_true")
output_group.add_argument('-o', '--out', help='prefix for output files (default: out.bayenv)', metavar="<STR>", type=str, action="store", default='bayenv')
args = parser.parse_args()

out_formats = ['b']
if (not args.bayenv and args.treemix) or args.bayenv:
	out_formats = ['b']

if args.treemix:
	out_formats.append('t')

if args.min_number:
	min_ind = args.min_number
else:
	min_ind = args.min_prop

#read in the vcf file
vcf = vcf.Reader(open(args.VCF, 'r'))


#create dictionary with population assignemtn for each individual
popdict = defaultdict(list)
samples = []
if args.popmap:
	if not os.path.isfile(args.popmap):
		print "The provided populationmap is not a valid file"
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
    			
			popdict[pop].append(ind)
			samples.append(ind)
	print "individuals per population:"
	for pop in sorted(popdict.keys()):
		print "%s:\t%s" %(pop, popdict[pop])
else:
	print "currently a populationmap is expected"
	sys.exit(0)

#check if the individuals in the populationmap fit the individuals in the vcf file
if not len(vcf.samples) == len(samples):
	print "vcf file does not contain the same number of samples (%i) specified in the populationmap" %len(samples)
	sys.exit(1)

for sam in vcf.samples:
	if not sam in samples:
		print "sample %s is not present in populationmap" %sam
		sys.exit(1)

###this produces a dictionary that contains the minimum number of indidividuals with data for each population
minimum_counts = {}
if isinstance(min_ind, float):
	for pop in popdict.keys():
		minimum_counts[pop] = int(len(popdict[pop])*min_ind)
else:
	for pop in popdict.keys():
		minimum_counts[pop] = min_ind

print "minimum counts of individuals with data:" 
for pop in sorted(minimum_counts.keys()):
	print "%s:\t%s" %(pop,minimum_counts[pop])
###############

#go through the actual data in the vcf file
popcounts = defaultdict(list)
SNPids = []
rem_coun = 0

for record in vcf: ## for each snp in the vcf
    SNPids.append(str(record.CHROM)+"\t"+str(record.POS)+"\t"+str(record.ID))
#    print record
#    print record.CHROM
#    print record.POS
#    print record.ID
    for pop in sorted(popdict.keys()):
#	print pop
	tempstring = ""
	for sample in popdict[pop]:
#		print sample
#		print str(record.genotype(sample)['GT'])

		tempstring = tempstring+str(record.genotype(sample)['GT'])

#	print tempstring		

	if tempstring.count("None") > 0:	#if there is at least one individual with missing data for the current locus
#		print "found missing data"
		if (len(popdict[pop]) - tempstring.count('None')) < minimum_counts[pop]:	#if there is less than the minimum number of individuals with data for this population
#			print "found missing data in %s individual(s) - break" % tempstring.count('None')
			minimum = len(popcounts[pop])	#record the number of loci currently recorded for the population that was identified to contain missing data for the current locus
#			print "population %s currrently contains: %s loci" %(pop,len(popcounts[pop]))
	
			check_list = []
			for pop in popcounts.keys():	#loop over all populations again
				if len(popcounts[pop]) > minimum:	#if there is a population that contains more loci than the the one that had been found to contain missing data for teh current locus. May happen if this population was succesfully processed for this locus before any missing data was encountert for a subsequent population
#					print "will pop last element from population %s" %pop
					popcounts[pop].pop()	#remove the last locus from the population
	
				check_list.append(len(popcounts[pop]))	
#				print "%s currrently contains: %s; %s" %(pop,popcounts[pop],len(popcounts[pop]))
	
#			print check_list
			SNPids.pop()
			rem_coun += 1
			break	#break out of the loop and go to next record

	refcount = str(tempstring.count('0'))
	altcount = str(tempstring.count('1'))
	popocount = refcount+","+altcount	
#	print popocount
	popcounts[pop].append(popocount)

#	print popcounts[pop][-1]

#for key in popcounts.keys():
#	print "%s (%s): %s" %(key, len(popcounts[key]), popcounts[key])

if (args.random_n == 0) or (len(popcounts[popcounts.keys()[0]]) <= args.random_n):
	rand = range(len(popcounts[popcounts.keys()[0]]))
#	outfile = "full_"+str(len(rand))+".bayenv"
	print "\nrequested full set of %i SNPs (%i loci were removed because they did not meet the minimum criteria)" % (len(rand), rem_coun)
else:
	rand = random.sample(range(len(popcounts[popcounts.keys()[0]])),args.random_n)
#	outfile = "random_"+str(len(rand))+".bayenv.in"
	print "\nrequested randomly selected subset of %i SNPs" % (len(rand))

##############

print "writing files:\n"
for f in out_formats:
	if f is 'b':
		print "\t%s.bayenv.SNPfile\n\t%s.bayenv.SNPmap\n" %(args.out, args.out)
		OUT = open(args.out+".bayenv.SNPfile","w")
		IDS = open(args.out+".bayenv.SNPmap","w")
		for i in rand:
			IDS.write(SNPids[i] + "\n")
			for j in range(2):
#				print j
				outlist = []
				outstring=""
				for pop in sorted(popdict.keys()):
					outlist.append(popcounts[pop][i].split(",")[j])
		
				outstring = "\t".join(outlist)
				outstring = outstring + "\t"
				OUT.write(outstring + "\n")

	elif f is 't':
		print "\t%s.treemix.SNPfile\n\t%s.treemix.SNPmap\n" %(args.out, args.out)
		OUT = open(args.out+".treemix.SNPfile","w")
		IDS = open(args.out+".treemix.SNPmap","w")
		outstring = " ".join(sorted(popdict.keys()))
		OUT.write(outstring + "\n")
		for i in rand:
			IDS.write(SNPids[i] + "\n")
			outlist = []
			for pop in sorted(popdict.keys()):
				outlist.append(popcounts[pop][i])

			outstring = " ".join(outlist)
			OUT.write(outstring + "\n")

	OUT.close()

sys.exit()
