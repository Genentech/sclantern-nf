#!/usr/bin/env python3

"""
This script takes the FLAG of each read in ORI_BAM and applies them to reads of the same QNAME in SNCR_BAM.

The input parameters are:
* ORI_BAM: Path of the original BAM file, BAM before splitting reads with GATK's SplitNCigarReads (SNCR) function.
* SNCR_BAM: Path of the split BAM file, BAM after splitting reads with SNCR.
* OUTPUT_BAM: Path of the output BAM file to be written (the '.bam' extension is added if missing).
* THREADS: Number of cores to use.

To run flagCorrection:
python3 flagCorrection.py \
  ORI_BAM \
  SNCR_BAM \
  OUTPUT_BAM \
  THREADS
"""

import sys
import os
import pysam
import multiprocessing
import argparse
import tempfile
import time
import subprocess

def main():
    start_time = time.time()

    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Apply flags from ORI_BAM to SNCR_BAM based on qname")
    parser.add_argument("ORI_BAM", help="Path to the original BAM file")
    parser.add_argument("SNCR_BAM", help="Path to the split BAM file")
    parser.add_argument("OUTPUT_BAM", help="Path to the output BAM file")
    parser.add_argument("THREADS", type=int, help="Number of cores to use")
    args = parser.parse_args()

    # Test input files and directories
    if not os.path.exists(args.ORI_BAM):
        sys.exit("File {} doesn't exist".format(args.ORI_BAM))
    if not os.path.exists(args.SNCR_BAM):
        sys.exit("File {} doesn't exist".format(args.SNCR_BAM))
    output_dir = os.path.dirname(os.path.abspath(args.OUTPUT_BAM))
    if not os.path.exists(output_dir):
        sys.exit("Directory {} doesn't exist".format(output_dir))

    # Create temporary directory
    temp_dir = tempfile.mkdtemp(prefix="flagCorrection_temp_dir_", dir=output_dir)

    # Load flags and qnames from the original BAM
    print("Loading flags and qnames from ORI_BAM")
    ori_bam = pysam.AlignmentFile(args.ORI_BAM, "rb")
    ori_qname_flag = {}
    for read in ori_bam.fetch(until_eof=True):
        ori_qname_flag[read.query_name] = read.flag
    ori_bam.close()

    # Get chromosome names from SNCR_BAM
    print("Retrieving chromosome names from SNCR_BAM")
    sncr_bam = pysam.AlignmentFile(args.SNCR_BAM, "rb")
    chr_names = sncr_bam.references
    sncr_bam.close()

    # Determine the number of threads to use
    threads = min(args.THREADS, len(chr_names))

    # Prepare tasks for parallel processing
    tasks = [(chr_name, args.SNCR_BAM, temp_dir, ori_qname_flag) for chr_name in chr_names]

    # Process chromosomes in parallel
    print("Processing chromosomes in parallel")
    pool = multiprocessing.Pool(processes=threads)
    temp_bam_files = pool.starmap(process_chromosome, tasks)
    pool.close()
    pool.join()

    # Ensure the output BAM filename ends with '.bam'
    if not args.OUTPUT_BAM.endswith('.bam'):
        args.OUTPUT_BAM += '.bam'

    # Merge temporary BAM files into the final output BAM
    print("Merging temporary BAM files into OUTPUT_BAM")
    samtools_command = ['samtools', 'merge', '-c', '-p', '-f', '-@', str(args.THREADS), args.OUTPUT_BAM] + temp_bam_files
    subprocess.check_call(samtools_command)

    # Clean up temporary files and directory
    for temp_file in temp_bam_files:
        os.remove(temp_file)
    os.rmdir(temp_dir)

    # Report time spent
    end_time = time.time()
    elapsed_time = (end_time - start_time) / 60  # in minutes
    print("\nflagCorrection finished after {:.2f} minutes.\n".format(elapsed_time))

def process_chromosome(chr_name, SNCR_BAM, temp_dir, ori_qname_flag):
    """
    Process a single chromosome: replace flags in SNCR_BAM with those from ORI_BAM based on qname.
    """
    # Open SNCR_BAM and create a temporary output BAM file for the chromosome
    sncr_bam = pysam.AlignmentFile(SNCR_BAM, "rb")
    temp_bam_filename = os.path.join(temp_dir, "{}_corrected.bam".format(chr_name))
    temp_bam = pysam.AlignmentFile(temp_bam_filename, "wb", template=sncr_bam)

    # Fetch reads for the chromosome and replace flags
    for read in sncr_bam.fetch(chr_name):
        qname = read.query_name
        if qname in ori_qname_flag:
            read.flag = ori_qname_flag[qname]
        else:
            # If qname not found in ori_bam, leave the flag unchanged
            pass
        temp_bam.write(read)

    sncr_bam.close()
    temp_bam.close()

    return temp_bam_filename

if __name__ == "__main__":
    main()