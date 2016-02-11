#!/usr/bin/python

import sys ## Exportable and executable
import os.path
from collections import defaultdict
import random
import argparse
from argparse import RawTextHelpFormatter
import vcf

parser = argparse.ArgumentParser(description='converts vcf format to SNPfile as required by Bayenv2 or Treemix')

parser.add_argument('VCF', help='vcf 4.0 (gzipped is supported)')
optional_group = parser.add_argument_group('Output options', 'The parameters in this group affect the format and name of the output')
parser.add_argument('-p', '--min_prop', help='minimum proportion of individuals with data per population (default: 1.00)', metavar="<FLOAT>", type=float, action="store", default=1.0)
parser.add_argument('-n', '--min_number', help='minimum number of individuals with data per population (default: 0, which means all individuals per population)', metavar="<INT>", type=int, action="store")
parser.add_argument('-r', '--random_n', help='number of randomly sampled loci to be outputted (default = 0, which means all loci)', metavar="INT", type=int, action="store", default=0)
parser.add_argument('-m', '--popmap', help='Tab delimited text file to assign individuals to populations. col1 = population ID, col2 = individual (as in vcf).', metavar="<FILE>", type=str, action="store")
parser.add_argument('-w', '--whitelist', help='Txt file containing a list of IDs for Loci to be retained.', metavar="<FILE>", type=str, action="store")
parser.add_argument('--pool', help='extract read counts as allele frequencies for each sample', action="store_true")
parser.add_argument('--min_global_MA_count', help='minimum number of global observations of minor allele', metavar="INT", type=int, action="store", default=0)
parser.add_argument('--min_global_MAF', help='minimum global minor allele frequency (MAF). This option will override --min_global_MA_count.', metavar="<FLOAT>", type=float, action="store", default=0)
parser.add_argument('--exclude_singletons', help='exclude any singleton loci, i.e. minor allele observed only once', action="store_true")
parser.add_argument('-v','--verbose', help='verbose', action="store_true")

output_group = parser.add_argument_group('Output options', 'The parameters in this group affect the format and name of the output')
output_group.add_argument('-b', '--bayenv', help='output bayenv format (default)', action="store_true")
output_group.add_argument('-t', '--treemix', help='output treemix format', action="store_true")
output_group.add_argument('-o', '--out', help='prefix for output files (default: out)', metavar="<STR>", type=str, action="store", default='out')
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

if args.exclude_singletons:
	args.min_global_MA_count=2

#read in the vcf file
vcf = vcf.Reader(open(args.VCF, 'r'))

###read in whitelist and check if it's ok
if args.whitelist:
	print "\nwhitelist specified\n"
	whitelist = defaultdict(int)
	WHITE = open(args.whitelist,"r")
	for wh in WHITE.readlines():
		whitelist[wh.strip()]+=1
	
	whitelist_count = len(whitelist)
	if not whitelist_count:
		print "\nwhitelist specified but file seems empty - might want to check that..\n"
		sys.exit(3)

#create dictionary with population assignemtn for each individual
popdict = defaultdict(list)
samples = []
if args.popmap:
	if not os.path.isfile(args.popmap):
		print "The provided populationmap is not a valid file"
		sys.exit(0)
	else:
		pops = open(args.popmap ,'r')
		if not args.pool:
			for pop in pops.readlines():
				pop = pop.strip()
    				ind = pop.split("\t")[0]
    				pop = pop.split("\t")[1]
				if not ind and pop:
					print "The populationmap is expected to have at least 2 tab delimited columns"
					sys.exit(0)
    				
				popdict[pop].append(ind)
				samples.append(ind)

			print "\nindividuals per population (as specified by populationmap):"
			for pop in sorted(popdict.keys()):
				print "%s (%i): %s" %(pop, len(popdict[pop]), ", ".join(popdict[pop]))

			#check if the individuals in the populationmap fit the individuals in the vcf file
			if not len(vcf.samples) >= len(samples):
				print "\nBe aware:\nvcf file contains %i samples, more than specified in the populationmap (%i)" %(len(vcf.samples), len(samples))
#				sys.exit(1)
	
			for sam in vcf.samples:
				if not sam in samples:
					print "sample %s is not present in populationmap" %sam
#					sys.exit(1)

			###this produces a dictionary that contains the minimum number of indidividuals with data for each population
			minimum_counts = {}
			if isinstance(min_ind, float):
				for pop in popdict.keys():
					minimum_counts[pop] = int(len(popdict[pop])*min_ind)
			else:
				for pop in popdict.keys():
					minimum_counts[pop] = min_ind

			print "\nspecified minimum counts of required individuals per population:" 
			for pop in sorted(minimum_counts.keys()):
				print "%s:\t%s" %(pop,minimum_counts[pop])
			###############

		elif args.pool:
			print "will interpret samples in the vcf file as populations"
			for pop in pops.readlines():
				pop = pop.strip()
    				
				popdict[pop].append(pop)
				samples.append(pop)

			print "populations (as specified in populationmap):"
			for pop in sorted(popdict.keys()):
				print "%s" %pop

else:
	print "please provide a populationmap"
	sys.exit(0)

if args.min_global_MAF:
	if args.min_global_MAF > 0.5:
		print "\nminimum MAF can't be higher than 0.5 - try again\n"
		sys.exit(4)
	print "\nyou specified minimum minor allele frequency (MAF): %f\n" %args.min_global_MAF
elif args.min_global_MA_count:
	print "\nyou specified minimum minor allele count: %i\n" %args.min_global_MA_count

#go through the actual data in the vcf file
popcounts = defaultdict(list)
SNPids = []
rem_coun = 0
list_count = 0

for record in vcf: ## for each snp in the vcf
	temp_pop = {}
	ref_count = 0 	#reference allele count per locus
	alt_count = 0	#alternative allele count per locus
#	print record.samples
	if args.whitelist:
#		print "currently %s elements in whitelist" %len(whitelist)
		if len(whitelist) == 0:
			break

		if not whitelist.has_key(str(record.ID)):
#			print "did not find %s in whitelist" %str(record.ID)
			continue	#if the current record is not present in the whitelist continue with the next iteration in the for loop, i.e. the next record
		else:
#			print "found %s in whitelist" %str(record.ID)
			del whitelist[str(record.ID)]

#	print record
#	print record.CHROM
#	print record.POS
#	print record.ID
	for pop in sorted(popdict.keys()):
#		print pop
		tempstring = ""
		if not args.pool:
			for sample in popdict[pop]:
#				print sample
#				print str(record.genotype(sample)['GT'])
	
				tempstring = tempstring+str(record.genotype(sample)['GT'])
	
#			print tempstring		
	
			if tempstring.count("None") > 0:	#if there is at least one individual with missing data for the current locus
#				print "found missing data"
				if (len(popdict[pop]) - tempstring.count('None')) < minimum_counts[pop]:	#if there is less than the minimum number of individuals with data for this population
#					print "found missing data in %s individual(s) - break" % tempstring.count('None')
##					minimum = len(popcounts[pop])	#record the number of loci currently recorded for the population that was identified to contain missing data for the current locus
#					print "population %s currrently contains: %s loci" %(pop,len(popcounts[pop]))
		
##					check_list = []
##					for pop in popcounts.keys():	#loop over all populations again
##						if len(popcounts[pop]) > minimum:	#if there is a population that contains more loci than the the one that had been found to contain missing data for teh current locus. May happen if this population was succesfully processed for this locus before any missing data was encountert for a subsequent population
#							print "will pop last element from population %s" %pop
##							popcounts[pop].pop()	#remove the last locus from the population
			
##						check_list.append(len(popcounts[pop]))	
#						print "%s currrently contains: %s; %s" %(pop,popcounts[pop],len(popcounts[pop]))
			
#					print check_list
##					SNPids.pop()
					rem_coun += 1
					break	#break out of the loop and go to next record
	
			refcount = str(tempstring.count('0'))
			altcount = str(tempstring.count('1'))
			popocount = refcount+","+altcount	
#			print popocount
#			popcounts[pop].append(popocount)
			temp_pop[pop] = popocount
			ref_count += int(refcount)
			alt_count += int(altcount)

#			print "refcount: %s\tref_count: %i" %(refcount,ref_count)
#			print "altcount: %s\talt_count: %i" %(altcount,alt_count)

		else:
#			for sample in popdict[pop]:
			if isinstance(record.genotype(pop)['AO'],list) or isinstance(record.genotype(pop)['AO'],list):
				rem_coun += 1
				break

			if not record.genotype(pop)['AO'] or not record.genotype(pop)['RO']:
#				print "population %s is missing data for locus %s, %s" %(pop, record.CHROM, record.POS)
#				print "The problematic population currently has %i elements" %len(popcounts[pop])
##				for pops in popcounts.keys():
#					print "pop: %s has %i elements" %(pops, len(popcounts[pops]))
##					if len(popcounts[pops]) > len(popcounts[pop]):
##						popcounts[pops].pop()
#						print "pop: %s has %i elements after removal of last element" %(pops, len(popcounts[pops]))
				rem_coun += 1
				break
			else:
#				print "%s: %s/%s" %(pop, record.genotype(pop)['RO'], record.genotype(pop)['AO'])
				ref_count += record.genotype(pop)['RO']
				alt_count += record.genotype(pop)['AO']
				temp_pop[pop] = str(record.genotype(pop)['RO'])+","+str(record.genotype(pop)['AO'])


	
	if len(temp_pop) == len(popdict):
		verbose_out="%s" %temp_pop
		if args.min_global_MAF:
			args.min_global_MA_count=int((ref_count+alt_count)*args.min_global_MAF)

		if not ref_count and not alt_count:
			verbose_out+="\texcluded: all loci are monomorphic"
			rem_coun += 1
#			continue
		elif (alt_count < args.min_global_MA_count):
			verbose_out+="\texcluded: minor allele only observed %s time(s) - minimum %i" %(alt_count, args.min_global_MA_count)
			rem_coun += 1
#			continue
		else:
			SNPids.append(str(record.CHROM)+"\t"+str(record.POS)+"\t"+str(record.ID))
			for key in temp_pop.keys():
				popcounts[key].append(temp_pop[key])
		if args.verbose:
			print "%s\t%s" %(str(record.CHROM)+"\t"+str(record.POS)+"\t"+str(record.ID),verbose_out)

if not args.whitelist:
	if (args.random_n == 0) or (len(popcounts[popcounts.keys()[0]]) <= args.random_n):
		rand = range(len(popcounts[popcounts.keys()[0]]))
#		outfile = "full_"+str(len(rand))+".bayenv"
		print "\nrequested full set of %i SNPs (%i loci were removed because they did not meet the minimum criteria)" % (len(rand), rem_coun)
	else:
		rand = random.sample(range(len(popcounts[popcounts.keys()[0]])),args.random_n)
#		outfile = "random_"+str(len(rand))+".bayenv.in"
		print "\nrequested randomly selected subset of %i SNPs" % (len(rand))
else:
	print "\nspecified whitelist with %i IDs" %whitelist_count
	print "\n%i (%.2f %%) of the loci specified in the whitelist were found" %(len(popcounts[popcounts.keys()[0]]), float(len(popcounts[popcounts.keys()[0]]))/whitelist_count*100)
	if len(popcounts[popcounts.keys()[0]]) != whitelist_count:
		print "The following loci were not present in the provided vcf file:"
		for key in sorted(whitelist.keys()):
			print key

	rand = range(len(popcounts[popcounts.keys()[0]]))
##############

print "\nwriting files:\n"
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
