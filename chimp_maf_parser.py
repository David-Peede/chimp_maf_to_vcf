#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 13:39:18 2021

@author: davidpeede
"""

import gzip
import sys


# sys.argv[1] = hg19.panTro6.synNet.maf.gz
# sys.argv[2] = chromosome to extract data for


def maf_reader(maf_file):
    """
    Reads in a .maf file and appends all lines in a synteny sequence alignment
    block to a list which will be parsed by the function maf_parser.
    """
    line = maf_file.readline()
    while line[0] != "a":
        line = maf_file.readline()
    block = []
    while line != "":
        line = maf_file.readline()
        if line == "" or line[0] == "a":
            yield block
            block = []
        elif line[0] == "s": block.append(line)

def maf_parser(maf_block):
    """
    Takes the output from the maf_reader function and creates a dictionary for
    each synteny sequence alignment block where the sequence name are the keys
    and the values are a dictionary consisting of the sequence start, size,
    strand oreintation, length, and the actual sequence.
    """
    output = {}
    for line in maf_block:
        source, start, size, strand, srcSize, seq = line.split()[1:]
        output[source] = {
                'start': int(start),
                'size': int(size),
                'strand': strand,
                'srcSize': int(srcSize),
                'seq': seq
                }
    return output

def extract_calls(maf_file,
                  ref_species='hg19',
                  contig='1',
                  min_length=1):
    """
    maf_file = bgziped .maf file (str)
    ref_species = reference source for the .maf file (str)
    contig = chromosome to extract allele calls for (str)
    min_length = minimum alignment length required to extract calls (int)
    ----------
    Outputs a .txt file to stdout in the format of:
    conting    pos    ref    target
    """
    chromosome = 'chr'+contig
    reference = ref_species+'.'+chromosome
    infile = gzip.open(maf_file, 'rt')
    outfile = sys.stdout
    allele = str.maketrans('-acgtn', 'NACGTN')
    blocks = maf_reader(infile)
    for block in blocks:
        block_info = maf_parser(block)
        sequences = block_info.keys()
        if reference != list(sequences)[0]:
            continue
        elif block_info[reference]['size'] <= min_length:
            continue
        else:
            true_len = block_info[reference]['size']
            align_len = len(block_info[reference]['seq'])
            ref_ind = [
                    i for i in range(align_len)
                    if block_info[reference]['seq'][i] != '-'
                    ]
            allele_calls = {}
            positions = range(
                    block_info[reference]['start']+1,
                    block_info[reference]['start']+1 + true_len,
                    )
            allele_calls[reference] = block_info[reference]['seq'].\
                replace('-', '').translate(allele)
            allele_calls[list(sequences)[1]] = block_info[list(sequences)[1]]['seq'].\
                translate(allele)
            for j in range(true_len):
                chrom_positions = '\t'.join([chromosome, str(positions[j])])
                refrence_calls = '\t'.join([allele_calls[list(sequences)[0]][j]])
                target_calls = '\t'.join([allele_calls[list(sequences)[1]][ref_ind[j]]])
                outfile.write(
                    chrom_positions\
                    + '\t'\
                    + refrence_calls\
                    + '\t'\
                    + target_calls\
                    + '\n')
    return

extract_calls(str(sys.argv[1]), 'hg19', str(sys.argv[2]))