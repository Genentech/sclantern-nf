#!/usr/bin/env python3

import argparse
import pysam
import numpy as np

def get_format(record,key):
    format_list = record.samples.items()[0][1].items()
    format_dict = dict(format_list)
    
    return(format_dict[key])

def TRGT_record_filter(record,motif_len = 1,depth = 10, al=1, ap = 0.5):
    coverage = get_format(record,"SD")
    if((coverage[0] is None) or ((coverage[0] < depth) or (coverage[1] < depth))):
        return(False)
        
    allele_length=get_format(record,"AL")
    if (allele_length[0] < al) or (allele_length[1] < al):
        return(False)
    allele_purity=get_format(record,"AP")
    if((allele_purity[0] is None) or ((allele_purity[0] < ap) or (allele_purity[1] < ap))):
        return(False)
        
    genotype = get_format(record,"GT")
    if((genotype[0] is None) or (genotype[0] == genotype[1])):
        return(False)
    
    motif = record.info["MOTIFS"]
    valid_motif_id = [len(i) >= motif_len for i in motif]
    if(sum(valid_motif_id) == 0):
        return(False)
    else:
        motif_repeat = get_format(record,"MC")
        motif_repeat = [i.split("_") for i in motif_repeat]
        flag = np.array(motif_repeat[0])[valid_motif_id] != np.array(motif_repeat[1])[valid_motif_id]
        if(sum(flag) == 0):
            return(False)
    return(True)
#TODO: ADD purity, length, and output motif
def TRGT_vcf_filter(vcf_file,out_file, motif_len = 1,depth=10, al=1, ap = 0.5):
    vcf = pysam.VariantFile(vcf_file,"r")
    out = pysam.VariantFile(out_file,"wb",header=vcf.header)

    out_list = []
    
    for record in vcf:
        if(TRGT_record_filter(record,motif_len,depth)):
            out_list.append(record)
            out.write(record)
    #pysam.tabix_index(out_file, preset="vcf", force=True)
    return(out_list)

def DV_record_filter(record,qual = 10,depth = 10):
    coverage = get_format(record,"DP")
    if((coverage is None) or (coverage < depth)):
        return(False)
    
    genotype = get_format(record,"GT")
    if((genotype[0] is None) or (genotype[0] == genotype[1])):
        return(False)
    
    if(record.qual < qual):
        return(False)
    
    return(True)

def DV_vcf_filter(vcf_file,out_file,qual = 10, depth = 10):
    vcf = pysam.VariantFile(vcf_file,"r")
    out = pysam.VariantFile(out_file,"wb",header=vcf.header)

    out_list = []
    
    for record in vcf:
        if(DV_record_filter(record,qual,depth)):
            out_list.append(record)
            out.write(record)
    #pysam.tabix_index(out_file, preset="vcf", force=True)
    return(out_list)

def main():
    parser = argparse.ArgumentParser(description="Filter TRGT VCF entries before variant calling")
    parser.add_argument("--in_vcf", required=True, help="Path to original VCF")
    parser.add_argument("--out_vcf", required=True, help="Path to filtered VCF")
    parser.add_argument("--mode", required=True, help="Either dv (for deepvariant) or trgt")
    args = parser.parse_args()

    if args.mode == "dv":
        DV_vcf_filter(args.in_vcf, args.out_vcf)
    elif args.mode == "trgt":
        TRGT_vcf_filter(args.in_vcf, args.out_vcf)
    else:
        raise Exception(f"Unrecognized mode {args.mode}: must be one of dv, trgt")

if __name__ == "__main__":
    main()
