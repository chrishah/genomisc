
# coding: utf-8

# In[ ]:

def normalize (csv, norm_prefix="", normalize=True, boxplot=False, boxplots_prefix=""):
    """
    The function parses a csv file and outputs normalized values
    and boxplots if desired
    """
    from collections import defaultdict
    import numpy as np
    import subprocess
    
    if normalize and not norm_prefix:
        raise IOError("You have to specify a prefix for the output files containing the normalized data - use 'norm_prefix='")
        
    if boxplot and not boxplots_prefix:
        raise IOError("You have to specify a prefix for the boxplot files to be produced - use 'boxplots_prefix='")
            
    populations = []    
    columns = {}
    indices = defaultdict(list)
    normalized = defaultdict(list)
    
    IDs = []
    pops = []
    
    INFILE = open(csv, 'r')

    headers = INFILE.readline().strip().split(",")
#    print headers
    headers.pop(0)
    IDs = sorted(headers)
    for env in headers:
        columns[env] = []
        
    for line in INFILE:
        line = line.strip()
        temp=line.split(",")
        populations.append(temp.pop(0))
#        print population
        for i in range(len(temp)):
#	   print temp[i]
#           print "%i:%s\n" %(i, header[i])
#	    if not temp[i] == 'NA':
	   try:
           	columns[headers[i]].append(float(temp[i]))
	   except ValueError:
            	columns[headers[i]].append(temp[i])
		
#        print columns
 
    pops = list(sorted(set(populations)))
#    print pops
    #find indexes for each
    for pop in pops:
#        print "finding indices for %s" %pop
        for i in range(len(populations)):
            if pop == populations[i]:
#                print i
                indices[pop].append(i)

#    print indices
    
#    print "\nCalculating means\n"
    for env in headers:
        per_pop = {}
	temp_per_env_list = []
	for value in columns[env]:
	     if not value == 'NA':
		temp_per_env_list.append(value)
#	print "global: %s\n" %(temp_per_env_list)
        for pop in pops:
#            print pop
#            print "%s - should be %i" %(env, len(indices[pop]))       
            per_pop_list = []
            for i in indices[pop]:
		if not columns[env][i] == 'NA':
                	per_pop_list.append(columns[env][i])
                
#            print per_pop_list
#            print "%s mean: %s" %(pop, np.mean(per_pop_list))
#            print "%s mean: %s" %(env, np.mean(columns[env]))
#            print "%s sd: %s" %(env, np.std(columns[env]))
            per_pop[pop] = per_pop_list

#            norm = (np.mean(per_pop_list) - np.mean(columns[env])) / np.std(columns[env])
#	    print "%s: %s\n" %(pop, per_pop_list)
            norm = (np.mean(per_pop_list) - np.mean(temp_per_env_list)) / np.std(temp_per_env_list)
#            print norm
            normalized[env].append(norm)
        if boxplot:
            print "Creating boxplot for %s\n" %env
            Rscript = boxplots_prefix+env+'.R'
            FH = open(Rscript, 'w')
            for pop in pops:
#                print per_pop[pop]
                vector = ""
                for v in per_pop[pop]:
                    vector+=repr(v)+','
                FH.write(pop+' <- c(%s)\n' %vector[:-1])
            FH.write("svg(filename = '"+boxplots_prefix+env+".svg')\n") #svg(filename = '$prefix-top-$cutoff.svg')
            FH.write("boxplot(%s, names = c('%s'), main = '%s')\n" %(", ".join(pops), "', '".join(pops), env)) 
            FH.write("dev.off()\n")
            FH.close()
            c = subprocess.Popen(['Rscript', Rscript], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            (output,err) = c.communicate()

            if err:
                print err
#            else:
#                print output
            
#    print normalized
    if normalize:
        print "\nnormalizing %s environmental factors across %s populations\nwriting to:\n\t%s.bayenv\n\t%s.csv" %(len(IDs), len(pops), norm_prefix, norm_prefix)
        OUTCSV = open(norm_prefix+'.csv', 'w')
        OUTCSV.write(",%s\n" %",".join(pops))
        OUTBAYENV = open(norm_prefix+'.bayenv',"w")
        for env in sorted(columns.keys()):
            outstring = ""
#            IDs.append(env)
            for n in normalized[env]:
                outstring += repr(n)+"\t"

    #        print outstring
            OUTBAYENV.write(outstring+"\n")
            OUTCSV.write("%s,%s\n" %(env, outstring.replace('\t',',')[:-1]))
        
        OUTBAYENV.close()
        OUTCSV.close()
    
    return pops, IDs
        


# In[ ]:

def split_for_Bayenv(infile, out_prefix):
    """
    This function takes a bayenv formatted multi-SNP file,
    splits it up into separate files (one SNP per file).
    """
    SNPcount = 0
    temp = []
    IN = open(infile, 'r')
    for line in IN:
        line = line.strip()
        temp.append(line)
        SNPcount+=1
        if (SNPcount % 2 == 0):
            OUT = open('%s-%07d.txt' %(out_prefix, SNPcount/2), 'w') #out_prefix+'-'+str(SNPcount/2)+'.txt', 'w')
#            print "%s\t\n" %"\n".join(temp)
            OUT.write("%s\t\n" %"\t\n".join(temp))
            OUT.close()
            temp = []
            
            


# In[ ]:

def calculate_rank_stats(SNP_map, infiles, ids, prefix): #these options are currently not implemented, threshold = 0.01, window = 50e3, sigma_factor = 3, bootstrap_rep = 100):
    """
    The function will calculate rank statistics (averages, standard deviations, etc)
    across a number of Bayenv replicates
    """
    #import 
    from collections import defaultdict
    import os
    import numpy as np
    
    #define global variables
    global_dict = defaultdict(dict)
    glob = defaultdict(dict) #holds all information for the SNPs
    reps = []
    extremes = {}
    rep_count = 0
    return_dict = {}
    
    #read in SNP_map
    SNPmap = open(SNP_map, 'r')
    i=0
    for line in SNPmap:
        global_dict[i]['chrom'] = line.strip()
        i += 1
    print "Total number of SNPs (according to the SNPmap): %i" %len(global_dict)
    
    #assess Bayenv files
    print "Number of Bayenv replicates: %i" %len(infiles)
    
    #assess environmental factors
    print "Number of environmental factors analysed: %i" %len(ids)

    #display settings
#    print "Sliding window settings:"
#    print "\tWindow size: %i bp" %window
#    print "\tSigma: %s" %sigma_factor

    #start processing
    print "parsing bayenv files"
    for bayenv_file in sorted(infiles):
        
        print "\nprocessing replicate %i:\n%s" %(rep_count, bayenv_file)
        reps.append(bayenv_file)
        fh = open(bayenv_file,'r')
        j=0 #this variable will hold the index of the current SNP
        sorting = defaultdict(dict) #This dictionary will contain the rank sorted SNPs
        for bf in fh:
            for factor_index in range(len(ids)):
#                print j
                factor = ids[factor_index]
#                print "processing factor %s" %factor
                if not glob.has_key(factor):
                    glob[factor] = defaultdict(dict)
            
                if not sorting[factor].has_key(float(bf.split('\t')[1+factor_index])):
                    sorting[factor][float(bf.split('\t')[1+factor_index])]=[j]
                else:
                    sorting[factor][float(bf.split('\t')[1+factor_index])].append(j)
            j += 1

        
        fh.close()
        
        for factor in sorting.keys(): #do sorting and add rank information to global dictionary
#            print factor
            rank = 0
            for r in sorted(sorting[factor].keys()):
##                print "rank: %i" %rank
##                print "bf: %f" %r
##                print sorting[r]
##                print len(sorting[r])
                for SNP in sorting[factor][r]:
##                    print SNP
                    if not glob[factor].has_key(SNP):
                        glob[factor][SNP]['ranks'] = [rank]
                    else:
                        glob[factor][SNP]['ranks'].append(rank)
                rank += len(sorting[factor][r])
        rep_count += 1
            
            
    output_columns = ['avg_rank', 'med_rank', 'p20_rank', 'std_rank', 'var_rank', 'mad_rank', 'avg_rank_rel', 'med_rank_rel', 'p20_rank_rel','var_rank_rel', 'var_rank_weight', 'var_weighted_avg_rank', 'var_weighted_rel_avg_rank']
    
    print "\nSUMMARY:\n"
    for fac in sorted(glob.keys()):
#        print "%s: %i" %(fac, len(glob[fac]))
        extremes[fac] = defaultdict(int)
        variances = []
#        for SNPid in glob[fac].keys()[0:10]:
#            print "%s: %s" %(SNPid, glob[fac][SNPid]['ranks'])
        for SNPid in glob[fac].keys():
#            print glob[fac][SNPid]['ranks']
#            print absolute_deviation_from_median(data=glob[fac][SNPid]['ranks'])
#            print find_max(absolute_deviation_from_median(data=glob[fac][SNPid]['ranks']))
            for maxi in find_max(absolute_deviation_from_median(data=glob[fac][SNPid]['ranks'])):
                extremes[fac][maxi] += 1
#            print extremes
            glob[fac][SNPid]['avg_rank'] = np.mean(glob[fac][SNPid]['ranks'])
            glob[fac][SNPid]['med_rank'] = np.median(glob[fac][SNPid]['ranks'])
            glob[fac][SNPid]['p20_rank'] = np.percentile(glob[fac][SNPid]['ranks'], 20)
            glob[fac][SNPid]['std_rank'] = np.std(glob[fac][SNPid]['ranks'])
            glob[fac][SNPid]['var_rank'] = np.var(glob[fac][SNPid]['ranks'])
            glob[fac][SNPid]['mad_rank'] = mad(data=glob[fac][SNPid]['ranks'])
            variances.append(glob[fac][SNPid]['var_rank'])
            glob[fac][SNPid]['avg_rank_rel'] = np.mean(glob[fac][SNPid]['ranks'])/len(glob[fac])
            glob[fac][SNPid]['med_rank_rel'] = np.median(glob[fac][SNPid]['ranks'])/len(glob[fac])
            glob[fac][SNPid]['p20_rank_rel'] = np.percentile(glob[fac][SNPid]['ranks'], 20)/len(glob[fac])
        #find maximum variance for the current factor
        max_var = max(variances)

        print "Wrting stats to %s" %(prefix+'_'+fac+'.txt')
        OUT = open(prefix+'_'+fac+'.txt','w')
        OUT.write("chrom\tbp\tSNPID\t"+"\t".join(output_columns)+'\n')

        for SNPid in glob[fac].keys(): #calculate relative and weighted variances
            glob[fac][SNPid]['var_rank_rel'] = glob[fac][SNPid]['var_rank'] / max_var
            glob[fac][SNPid]['var_rank_weight'] = 1-glob[fac][SNPid]['var_rank_rel']
            glob[fac][SNPid]['var_weighted_avg_rank'] = glob[fac][SNPid]['var_rank_weight'] * glob[fac][SNPid]['avg_rank']
            glob[fac][SNPid]['var_weighted_rel_avg_rank'] = glob[fac][SNPid]['var_rank_weight'] * glob[fac][SNPid]['avg_rank_rel']

            temp_list = []
#            outstring = str(SNPid)+','
            
            outstring = global_dict[SNPid]['chrom']+'\t'
            for column in output_columns:
#                print column
#                print glob[fac][SNPid][column]
                temp_list.append(str(glob[fac][SNPid][column]))
            outstring += "\t".join(temp_list)
            OUT.write(outstring+'\n')

        OUT.close()
        

        counts = [extremes[fac][i] for i in sorted(extremes[fac])]
        perc = [float(counts[i])/len(glob[fac])*100 for i in range(len(counts))]
        ex = find_max([extremes[fac][i] for i in sorted(extremes[fac])])
        print "factor %s\treplicate %s gave the most extreme ranks for %.2f %% of the SNPs" %(fac, ex[0], perc[ex[0]])

    return_dict['global'] = glob
    return_dict['extremes'] = extremes
    return_dict['SNPids'] = global_dict
    return return_dict

# In[ ]:

def mad(data):
    """
    find the 'median absolute deviation'
    """
    import numpy as np
   
    return np.median(np.abs(data - np.median(data)))


# In[ ]:

def absolute_deviation_from_median(data):
    """
    find the absolute deviations from median for a list of values
    """
    import numpy as np
    
    mads = []
    med = np.median(data)
    for d in data:
        mads.append(np.abs(d-med))
        
    return mads


# In[ ]:

def find_max(data):
    """
    find the index of the maximum value in a list
    """
    from collections import defaultdict

    d = defaultdict(list)
    for i, x in enumerate(data):
        d[x].append(i)

    k = max(d.keys())
    return d[k]


# In[ ]:

def plot_pope(files_list, cutoff, num_replicates):
    """
    The function calls a shell script that configures and runs an 
    R script to plot pope plots and extract the top SNPs
    """
    import subprocess
    
    for ID in sorted(files_list):
        print "processing %s:" %ID[0]
        print "data in file: %s" %ID[1]
        c = subprocess.Popen(['sh','plot_pope.sh',ID[1],str(cutoff),ID[0],str(num_replicates)], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        (output,err) = c.communicate()

        if err:
            print err
        else:
            print output
        


# __Define the function__ that excludes the most extreme replicate.

# In[ ]:

def exclude_extreme_rep(dictionary, ids, prefix, cutoff=0):
    """
    Calculate rank statistics while excluding the one most extreme replicate for every factor
    """
    from collections import defaultdict
    import numpy as np
    
    output_columns = ['avg_rank', 'med_rank', 'std_rank', 'var_rank', 'mad_rank', 'avg_rank_rel', 'med_rank_rel', 'var_rank_rel', 'var_rank_weight', 'var_weighted_avg_rank', 'var_weighted_rel_avg_rank']
    glob = defaultdict(dict)
    
    for fac in sorted(ids):
        counts = [dictionary['extremes'][fac][i] for i in sorted(dictionary['extremes'][fac])]
        perc = [float(counts[i])/len(dictionary['global'][fac])*100 for i in range(len(counts))]
        ex = find_max([dictionary['extremes'][fac][i] for i in sorted(dictionary['extremes'][fac])])
        print "\nfactor %s\treplicate %s gave the most extreme ranks for %.2f %% of the SNPs" %(fac, ex[0], perc[ex[0]])

        if perc[ex[0]] > cutoff*100:
            print "will re-calculate stats without replicate %s" %ex[0]
        else:
            print "none of the replicates exceeds the %s threshold" %cutoff
            continue
        
        glob[fac] = defaultdict(dict)
        variances = []
        for SNPid in dictionary['global'][fac].keys():
#            print dictionary['global'][fac][SNPid]['ranks']
            temp_list = dictionary['global'][fac][SNPid]['ranks'][:]
            temp_list.pop(ex[0])
#            print temp_list
            glob[fac][SNPid]['avg_rank'] = np.mean(temp_list)
            glob[fac][SNPid]['med_rank'] = np.median(temp_list)
            glob[fac][SNPid]['std_rank'] = np.std(temp_list)
            glob[fac][SNPid]['var_rank'] = np.var(temp_list)
            glob[fac][SNPid]['mad_rank'] = mad(data=temp_list)
            variances.append(glob[fac][SNPid]['var_rank'])
            glob[fac][SNPid]['avg_rank_rel'] = np.mean(temp_list)/len(dictionary['global'][fac])
            glob[fac][SNPid]['med_rank_rel'] = np.median(temp_list)/len(dictionary['global'][fac])
        #find maximum variance for the current factor
        max_var = max(variances)

        print "Wrting stats to %s" %(prefix+'_'+fac+'_ex_rep_'+str(ex[0])+'.txt')
        OUT = open(prefix+'_'+fac+'_ex_rep_'+str(ex[0])+'.txt','w')
        OUT.write("chrom\tbp\tSNPID\t"+"\t".join(output_columns)+'\n')

        for SNPid in glob[fac].keys(): #calculate relative and weighted variances
            glob[fac][SNPid]['var_rank_rel'] = glob[fac][SNPid]['var_rank'] / max_var
            glob[fac][SNPid]['var_rank_weight'] = 1-glob[fac][SNPid]['var_rank_rel']
            glob[fac][SNPid]['var_weighted_avg_rank'] = glob[fac][SNPid]['var_rank_weight'] * glob[fac][SNPid]['avg_rank']
            glob[fac][SNPid]['var_weighted_rel_avg_rank'] = glob[fac][SNPid]['var_rank_weight'] * glob[fac][SNPid]['avg_rank_rel']

            temp_list = []
#            outstring = str(SNPid)+','
            outstring = dictionary['SNPids'][SNPid]['chrom']+'\t'
            for column in output_columns:
#                print column
#                print glob[fac][SNPid][column]
                temp_list.append(str(glob[fac][SNPid][column]))
            outstring += "\t".join(temp_list)
            OUT.write(outstring+'\n')

        OUT.close()


# In[ ]:

def parse_gff(gff):
    """
    parse gff file
    """
    
    gff_dict = {}
    
    gff_fh = open(gff,'r')
    for line in [line.strip() for line in gff_fh]:
#        print line.split('\t')
        if line.split('\t')[2] == 'CDS':
            gene = line.split('\t')[8].split(' ')[3].replace('"','').replace(';','') #This line needs to be adujsted to the gff format
            if not gff_dict.has_key(line.split('\t')[0]):
                gff_dict[line.split('\t')[0]] = {}
            gff_dict[line.split('\t')[0]][line.split('\t')[3]] = gene
            gff_dict[line.split('\t')[0]][line.split('\t')[4]] = gene
            
    return gff_dict    
        


# In[ ]:

def find_genes(rank_stats, gff, distance):
    """
    find genes up and downstream of SNPs
    """    
    from collections import defaultdict
    candidates = defaultdict(dict)
    return_dict = defaultdict(dict)
    
    for tsv in rank_stats:
        print "processing rank statistic file: %s" %tsv
        ID = tsv.split('/')[-1].replace('.tsv','')
        candidates[ID] = defaultdict(list)
        rank_stats_fh = open(tsv, 'r')
        rank_stats_fh.readline()
        for SNP in [SNP.strip() for SNP in rank_stats_fh]:
            rank_elem = SNP.split('\t')
            if not candidates[ID].has_key(rank_elem[0]):
                candidates[ID][rank_elem[0]] = []
            candidates[ID][rank_elem[0]].append(rank_elem[1:3])
#            print candidates[tsv][rank_elem[0]]
            

    for tsv in sorted(candidates.keys()):
        print "%s:" %tsv
        for chrom in sorted(candidates[tsv]):
#            print chrom
            for hot in candidates[tsv][chrom]:
                gene_list = []
                temp = []
                nr_genes = []
                lower = int(hot[0])-(distance*1000)
                upper = int(hot[0])+(distance*1000)
#                print "looking at %s" %hot[0]
                if not gff.has_key(chrom):
#                    print "no genes found on %s\n" %chrom
                    continue                
                else:
                    for pos in gff[chrom].keys():
                        temp.append(int(pos))
                    
                for pos in sorted(temp):
                    if pos >= lower and pos <= upper:
#                        print pos,gff[chrom][str(pos)]
                        gene_list.append(gff[chrom][str(pos)])
                    elif pos > upper:
                        break
                 
                nr_genes = list(set(gene_list))
                for unique_gene in nr_genes:
#                    print [chrom,hot[0],hot[1],unique_gene]
                    if not return_dict[tsv].has_key('genes'):
                        return_dict[tsv]['columns'] = ['chrom','bp','ID','gene']
                        return_dict[tsv]['genes'] = []
                    return_dict[tsv]['genes'].append([chrom,hot[0],hot[1],unique_gene])

        if not return_dict.has_key(tsv):
            return_dict[tsv]['genes'] = []
            return_dict[tsv]['columns'] = ['chrom','bp','ID','gene']
        print "identified %i gene(s)" %len(return_dict[tsv]['genes'])
        
    return return_dict


# In[ ]:

def annotate_genes(SNPs_to_genes, annotations, whitelist=[]):
    """
    fetch annotation for genes from file produced by Blast2GO
    """
    from collections import defaultdict
    
    annotation = defaultdict(list)
    
    if whitelist:
        for id in whitelist:
            if not SNPs_to_genes.has_key(id):
                raise IOError("You provide an analysis id %s that is not in the dictionary" %id)        
    else:
        whitelist = SNPs_to_genes.keys()[:]

    anno_fh = open(annotations, 'r')
    header = anno_fh.readline().strip().split('\t')
    annotation['header'] = header[1:]
    annotation['genes'] = defaultdict(list)
    for line in [line.strip() for line in anno_fh]:
        annotation['genes'][line.split('\t')[0]] = line.split('\t')[1:]
        
#    for gene in annotation['genes'].keys()[:10]:
#        print gene,annotation['genes'][gene]

    for analysis_id in whitelist:
        print analysis_id
        if len(SNPs_to_genes[analysis_id]['genes']) > 0:
            print "adding annoation for %s" %analysis_id
            for index in range(len(SNPs_to_genes[analysis_id]['genes'])):
                if annotation['genes'].has_key(SNPs_to_genes[analysis_id]['genes'][index][-1]):
                    
                    if len(SNPs_to_genes[analysis_id]['columns']) == 4:
#                        print SNPs_to_genes[analysis_id]['columns']
#                        print "extend the headers"
                        SNPs_to_genes[analysis_id]['columns'].extend(annotation['header'])
#                        print SNPs_to_genes[analysis_id]['columns']
#                    print annotation['genes'][SNPs_to_genes[analysis_id]['genes'][index][-1]]
                    SNPs_to_genes[analysis_id]['genes'][index].extend(annotation['genes'][SNPs_to_genes[analysis_id]['genes'][index][-1]])
                elif len(SNPs_to_genes[analysis_id]['genes'][index]) == 4 and not annotation['genes'].has_key(SNPs_to_genes[analysis_id]['genes'][index][-1]):
                    print "no annoation found for %s" %SNPs_to_genes[analysis_id]['genes'][index][-1]
            
    
        else:
            print "nothing to annotate - 0 candidate genes identified for %s" %analysis_id
    


# In[ ]:

def write_candidates(SNPs_to_genes, whitelist=[], rename=[], out_dir='./'):
    """
    write out SNP to candidate genes text files (will be named *.genes.tsv)
    """
    
    if rename:
        if not len(rename) == len(whitelist):
            raise IOError("If you provide a list with new names it needs to be the same length as the whitelist")

    if whitelist:
        for id in whitelist:
            if not SNPs_to_genes.has_key(id):
                raise IOError("You provide an analysis id %s that is not in the dictionary" %id)        
    else:
        whitelist = SNPs_to_genes.keys()[:]
                
    for id in sorted(whitelist):
        print id
        if not len(SNPs_to_genes[id]['genes']) >= 1:
            print "0 candidate genes found"
            continue
        else:
            if len(SNPs_to_genes[id]['columns']) == 4:
                print "writing to: %s" %(out_dir+id+'.genes.tsv')
                out_fh = open(out_dir+id+'.genes.tsv','w')
            else:
                print "writing to: %s" %(out_dir+id+'.genes.annotated.tsv')
                out_fh = open(out_dir+id+'.genes.annotated.tsv','w')

            out_fh.write("%s\n" %"\t".join(SNPs_to_genes[id]['columns']))
            for gene in SNPs_to_genes[id]['genes']:
                out_fh.write("%s\n" %"\t".join(gene))
            
        out_fh.close()


# In[ ]:



