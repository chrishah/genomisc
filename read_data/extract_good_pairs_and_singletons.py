#!/usr/bin/python
"""This script takes a fastq file and splits good pairs and orphans
Usage:
	extract_good_pairs_and_singletons.py fastq_file    
"""


import sys

def open_fastq(filename, r=False, w=False):
    """
    Function that opens a fastq file handle
    """
    
    import gzip
    
    if not r and not w:
        raise IOError('Either reading (r) or writing (w) needs to be specified\n')
        
    if r and w:
        raise IOError('Choose either reading (r) or writing (w) \n')
    
    if r:
        mode = 'r'
    elif w:
        mode = 'w'
    
    if filename.endswith('.gz'):
        FH = gzip.open(filename, mode+'b')
    else:
        FH = open(filename, mode)
        
    return FH

def extract_good_pairs_and_singletons(to_process, out_dir):
    """
    The function parses a fastq file (gzipped supported)
    and separates paired end from singleton reads
    """
    import gzip
    from Bio import SeqIO
    from collections import defaultdict
    
    id_dict = defaultdict(list)
    
    to_process_FH = open_fastq(to_process, r=True)
    pe_1_FH = open_fastq(out_dir+'/pe-1.fastq.gz', w=True)
    pe_2_FH = open_fastq(out_dir+'/pe-2.fastq.gz', w=True)
    se_FH = open_fastq(out_dir+'/se.fastq.gz', w=True)
    
    for read in SeqIO.parse(to_process_FH, 'fastq'):
        ID = read.id[:-2]
        id_dict[ID].append(read)
        if len(id_dict[ID]) == 2:
            id_dict[ID] = sorted(id_dict[ID])
            SeqIO.write(id_dict[ID][0], pe_1_FH, 'fastq')
            SeqIO.write(id_dict[ID][1], pe_2_FH, 'fastq')
            del id_dict[ID]
            
    for ID in id_dict.keys():
        SeqIO.write(id_dict[ID][0], se_FH, 'fastq')
               
    to_process_FH.close()
    pe_1_FH.close()
    pe_2_FH.close()
    se_FH.close()
    del id_dict

if __name__ == '__main__':
    try:
        file1 = sys.argv[1]
    except:
        print __doc__
        sys.exit(1)

    extract_good_pairs_and_singletons(to_process=sys.argv[1], out_dir='./')
