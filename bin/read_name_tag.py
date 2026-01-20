import pysam
import argparse
import os

def add_tags(bam_file,out_file,tags,suffix = "",sep = "_"):
    bam = pysam.AlignmentFile(bam_file, "rb",check_sq=False)
    header = bam.header.copy()
    
    out = pysam.AlignmentFile(out_file, "wb", header=header)
    
    for read in bam:
        qname = read.qname
        flag = True
        for i in tags:
            if read.has_tag(i):
                qname = qname + sep+read.get_tag(i)
            else:
                flag = False
        qname = qname+sep+suffix
        read.qname = qname
        if flag:
            out.write(read)
    
    bam.close()
    out.close()
    
    return(out_file)
