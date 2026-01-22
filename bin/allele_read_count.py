#!/usr/bin/env python

import argparse

import re

import pysam
from Bio import Align
import numpy as np
import pandas as pd

def sub_reads(read,start,end):
    pos = read.pos
    if(pos > start):
        return(None)
    vec = []
    for i in read.cigartuples:
        id = []
        if(i[0] in [0,7,8]):
            id = list(range(pos,pos+i[1]))
            pos = pos+i[1]-1
        elif(i[0] == 1):
            id = list(np.repeat(pos,i[1]))
            pos = pos+1
        elif(i[0]==2):
            pos = pos+i[1]+1
        vec = vec+id
        if pos > end:
            break
    vec = np.array(vec)
    if(np.max(vec) < end):
        return(None)
    indices = np.where((vec >= start) & (vec <= end))[0]
    
    out = ""
    if(len(indices) > 0):
        out = read.query_alignment_sequence[min(indices):max(indices)]
    
    return(out)

def ref_to_alt(seq,ref,alts,pos):
    out = [seq]
    if alts is None:
        alts = [""]
    for i in alts:
        out = out+[seq[0:pos]+i+seq[(pos+len(ref)):]]
    out = sorted(out, key=len, reverse=True)
    return(out)

def aligner_init(match = 2,mismatch = -1,internal_gap = -1, end_gap = -0.5,mode = "global"):
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.match_score = match
    aligner.mismatch_score = mismatch
    aligner.internal_gap_score = internal_gap
    aligner.end_gap_score = end_gap
    # aligner.open_gap_score  = open_gap
    # aligner.extend_gap_score = extend_gap
    return(aligner)


def loci_assign_reads(bam, fasta, chrom,pos,ref,alts,flank = 0):
    alted_ref = fasta.fetch(start = pos-flank-1,end = pos+len(ref)+flank-1,region = chrom)
    alted_ref = ref_to_alt(alted_ref,ref,alts,flank)

    # result = pd.DataFrame({'cell': [], 'UMI': [], "allele":[]})
    result = pd.DataFrame({'qname': [], "allele":[],"score":[]})
    aligner = aligner_init()
    for read in bam.fetch(chrom, pos-flank-1, pos+len(ref)+flank-1):
        subseq = sub_reads(read,pos-flank-1, pos+len(ref)+flank-1)
        if(subseq == None):
            continue
#         print(subseq)
        scores = []
        for i in alted_ref:
            temp = -1
            if((len(subseq) > 0) & (len(i) > 0)):
                alignments = aligner.align(subseq,i)
                temp = alignments.score
            else:
                temp = -abs(len(i)-len(subseq))
            scores.append(temp)
#         print(scores)
#         print(alted_ref)
        if(sum(np.array(scores) > 0) > 0):
            id = np.where(np.array(scores)==max(scores))[0][0]
            # new_row = {'cell': read.get_tag("CB"), 'UMI': read.get_tag("XM"),"allele":alts[id]}
            new_row = {'qname': read.qname,"allele":alted_ref[id],"score":max(scores)}
            result = result._append(new_row, ignore_index=True)
    # result = result.groupby(["cell","UMI","allele"]).size().reset_index(name='count')
    result["chrom"] = chrom
    result["pos"] = pos 
    return(result)


#def allele_count(
def assign_reads(
        bam_file, genome_fasta, vcf_file, flank = 0,
):
    bam = pysam.AlignmentFile(bam_file, "rb")
    genome= pysam.FastaFile(genome_fasta)
    vcf = pysam.VariantFile(vcf_file)
    
    result = []

    for record in vcf:
        sub = loci_assign_reads(bam,genome,record.chrom,record.pos,record.ref,record.alts,flank)
        result.append(sub)
    result = pd.concat(result, ignore_index=True)

    return result

# df_assignments should be a dataframe returned by assign_reads()
def assignments_to_counts(
        df_assignments,
        re_qname=r"molecule/\d+_([ACTG]+)_([ACTG]+)_.*",
        cb_group=1, umi_group=2
):
    cb = []
    umi = []
    loc = []

    re_qname = re.compile(re_qname)
    for _, row in df_assignments.iterrows():
        matched = re_qname.match(row['qname'])
        cb.append(matched.group(cb_group))
        umi.append(matched.group(umi_group))

        loc.append(f"{row['chrom']}_{row['pos']}")

    df_assignments['cell'] = cb
    df_assignments['umi'] = umi
    df_assignments['locus'] = loc

    return(
        df_assignments.groupby(['cell', 'umi', 'allele', 'locus'])
        .size()
        .reset_index(name='count')
    )

def main():
    parser = argparse.ArgumentParser(description="Call variants on UMIs and summarize read counts.")
    parser.add_argument("--bam_file", required=True, help="Input bam file.")
    parser.add_argument("--genome_fasta", required=True, help="Reference fasta.")
    parser.add_argument("--variants_vcf", required=True, help="VCF used for variant positions and alleles.")
    parser.add_argument("--out_count_path", required=True, help="Output csv path.")
    parser.add_argument("--sample_name", required=False, default="", help="Add sample name as prefix to cell barcode")
    args = parser.parse_args()

    df_assignments = assign_reads(
        args.bam_file,
        args.genome_fasta,
        args.variants_vcf,
    )

    counts = assignments_to_counts(df_assignments)

    if args.sample_name:
        counts['cell'] = args.sample_name + "_" + counts['cell']

    counts.to_csv(
        args.out_count_path,
        sep='\t', index=False
    )

if __name__ == "__main__":
    main()
