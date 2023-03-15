#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 19 13:39:18 2021
@author: davidpeede
"""

import gzip
import sys

def maf_output_to_VCF(maf_parser_output, header_file):
    """
    maf_parser_output = bgzipped output from the extract_calls function from
    updated_maf_parser.py
    header_file = a txt file with a generic VCF header.
    ----------
    Outputs a .txt file to stdout in the VCF format
    """
    outfile = sys.stdout
    with open(header_file) as infile:
        for line in infile:
            outfile.write(line)
    outfile.write('\n')
    with gzip.open(maf_parser_output, 'rt') as infile:
        for line in infile:
            spline = line.split()
            if spline[2] == 'N' or spline[3] == 'N':
                continue
            elif spline[3] == spline[2]:
                geno = '0|0'
                ALT = '.'
            else:
                geno = '1|1'
                ALT = spline[3]
            CHR = spline[0][3:]
            POS = spline[1]
            ID = '.'
            REF = spline[2]
            QUAL = '.'
            FILTER = 'PASS'
            INFO = '.'
            FORMAT = 'GT'
            outfile.write(
                    CHR\
                    + '\t'\
                    + POS\
                    + '\t'\
                    + ID\
                    + '\t'\
                    + REF\
                    + '\t'\
                    + ALT\
                    + '\t'\
                    + QUAL\
                    + '\t'\
                    + FILTER\
                    + '\t'\
                    + INFO\
                    + '\t'\
                    + FORMAT\
                    + '\t'\
                    + geno\
                    + '\n')
    return

maf_output_to_VCF(str(sys.argv[1]), './generic_vcf_header.txt')