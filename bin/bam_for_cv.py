#!/usr/bin/env python

import pysam
from Bio import Align
import numpy as np
import pandas as pd
import argparse
import warnings
from multiprocessing import Pool
import os
import uuid
import pathlib
import subprocess

from read_name_tag import add_tags

def args_parse():
    parser = argparse.ArgumentParser(description='Prepare the scRNA-seq bam file for variant calling')
    parser.add_argument('--sample_name', required=True, help="Sample name to be appended to read qnames")
    parser.add_argument('--bam', '-b', required=True, help='Path to the input bam file')
    parser.add_argument('--outbam', '-o', required=True, help='Path to tagged output bam')
    parser.add_argument('--minimap2', '-m', default="minimap2", help='Path to the minimap2 tool')
    parser.add_argument('--samtools', '-s', default="samtools", help='Path to the samtools tool')
    parser.add_argument('--genome', '-g', required=True, help='Path to the reference genome fasta')
    parser.add_argument('--CB', '-c', default="CB", help='The tag name of the cell barcode in the bam file')
    parser.add_argument('--UMI', '-u', default="XM", help='The tag name of the UMI in the bam file')
    parser.add_argument('--thread', '-t', type=int, default=1, help='The number of threads to use')

    return parser.parse_args()


def bam_flag(bam_file,CB_tag = "CB",UMI_tag = "XM",ncheck = 10):
    bam = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
    
    CB_flag = False
    UMI_flag = False
    qual_flag = False
    
    i = 0
    for read in bam:
        if(read.has_tag(CB_tag)):
            CB_flag = True
        if(read.has_tag(UMI_tag)):
            UMI_flag = True
        if(read.qual != None):
            qual_flag = True
        if CB_flag and UMI_flag and qual_flag:
            break
        i = i + 1
        if(i > ncheck):
            break
    bam.close()
    return([CB_flag,UMI_flag,qual_flag])

def bam_process(sample_name, bam_file,
                out_path, genome,
                minimap2, samtools,
                CB_tag = "CB",UMI_tag = "XM",nthread = 1):
    flag = bam_flag(bam_file)
    
    if(not flag[0]):
        message = "The cell barcode hasn't been identified for each read, please do this first!"
        RuntimeError(message)
    if(not flag[0]):
        message = "The UMI tag is not found for each read, please check the tag name, otherwise each read would be assigned with a distinct UMI."
        warnings.warn(message)
    
    tags = np.array([CB_tag,UMI_tag])[np.array(flag[0:2])]
    add_tags(bam_file, out_path, tags, suffix = sample_name)
    
    if(not flag[2]):
        message = "There is no sequence quality, would do remapping to generate a uniform sequence quality"
        warnings.warn(message)
    
args = args_parse()
print(args)

if __name__ == '__main__':
    args = args_parse()

    bam_process(
        args.sample_name,
        args.bam, args.outbam, args.genome, 
        args.minimap2, args.samtools, 
        args.CB, args.UMI, args.thread
    )
