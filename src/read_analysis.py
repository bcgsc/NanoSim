#!/usr/bin/env python

"""
Created on Apr 10, 2015

@author: Chen Yang

This script generates read profiles Oxford Nanopore 2D reads.

"""


from __future__ import print_function
from __future__ import with_statement
from subprocess import call
from time import strftime
import sys
import os
import getopt
import numpy
import head_align_tail_dist as align
import get_besthit
import besthit_to_histogram as error_model


# Usage information
def usage():
    usage_message = "./read_analysis.py <options>\n" \
                    "<options>: \n" \
                    "-h : print usage message\n" \
                    "-i : training ONT real reads, must be fasta files\n" \
                    "-r : reference genome of the training reads\n" \
                    "-m : User can provide their own alignment file, with maf extension\n" \
                    "-o : The prefix of output file, default = 'training'\n"

    sys.stderr.write(usage_message)


def main(argv):
    # Parse input and output files
    infile = ''
    outfile = 'training'
    ref = ''
    maf_file = ''
    model_fit = True
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:m:", ["infile=", "ref=", "outfile=", "model_fit="])
    except getopt.GetoptError:
        usage()
        sys.exit(1)

    for opt, arg in opts:
        if opt == '-h':
            usage()
            sys.exit(0)
        elif opt in ("-i", "--infile"):
            infile = arg
        elif opt in ("-r", "--ref"):
            ref = arg
        elif opt == "-m":
            maf_file = arg
        elif opt in ("-o", "--outfile"):
            outfile = arg
        elif opt == "-model_fit":
            model_fit = arg
        else:
            usage()
            sys.exit(1)

    if infile == '' or ref == '':
        print("Please specify the training reads and its reference genome!")
        usage()
        sys.exit(1)

    # Read pre-process
    in_fasta = outfile + ".fasta"
    if in_fasta == infile:
        in_fasta = outfile + "_processed.fasta"
    out_fasta = open(in_fasta, 'w')
    dic_reads = {}
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '>':
                name = '-'.join(line.strip()[1:].split())
                dic_reads[name] = ""
            else:
                dic_reads[name] += line.strip()
    for k, v in dic_reads.items():
        out_fasta.write('>' + k + '\n' + v + '\n')
    out_fasta.close()

    # ALIGNMENT if maf file is not provided
    out_maf = outfile + ".maf"
    if out_maf == maf_file:
        out_maf = outfile + "_processed.maf"
    if maf_file != '':
        call("grep '^s ' " + maf_file + " > " + out_maf, shell=True)
    else:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
        call("lastdb ref_genome " + ref, shell=True)
        call("lastal -a 1 ref_genome " + in_fasta + " | grep '^s ' > " + out_maf, shell=True)

    # LENGTH DISTRIBUTION
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": length distribution analysis\n")
    # Get besthit of alignment
    unaligned_length = get_besthit.besthit_and_unaligned(in_fasta, out_maf, outfile)

    # Generate the histogram of aligned reads
    num_aligned = align.head_align_tail(outfile)

    # Length distribution of unaligned reads
    out1 = open(outfile + "_unaligned_length_ecdf", 'w')
    out2 = open(outfile + "_unaligned_length.txt", 'w')
    out2.write('\n'.join(str(x) for x in unaligned_length))
    out2.close()

    num_unaligned = len(unaligned_length)
    if num_unaligned != 0:
        max_length = max(unaligned_length)
        hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length, 50), density=True)
        cdf = numpy.cumsum(hist_unaligned * 50)
        out1.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
        out1.write("bin\t0-" + str(max_length) + '\n')
        for i in xrange(len(cdf)):
            out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    else:
        out1.write("Aligned / Unaligned ratio:\t100%\n")

    out1.close()

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile)

    if model_fit:
        path = sys.argv[0].split("/")
        R_path = '/'.join(path[:-1]) + '/' + "model_fitting.R"
    	if os.path.isfile(R_path):
            call("R CMD BATCH '--args prefix=\"" + outfile + "\"' " + R_path, shell=True)
        else:
            print("Could not find 'model_fitting.R' in ../src/\nMake sure you copied the whole source files from Github.")

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main(sys.argv[1:])

