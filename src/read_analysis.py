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
try:
    from six.moves import xrange
except ImportError:
    pass
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
                    "-b : number of bins (for development), default = 20\n" \
                    "-o : The prefix of output file, default = 'training'\n" \
                    "--no_model_fit : Skip the model fitting step\n"

    sys.stderr.write(usage_message)


def main(argv):
    # Parse input and output files
    infile = ''
    outfile = 'training'
    ref = ''
    maf_file = ''
    model_fit = True
    num_bins = 20
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:m:b:", ["infile=", "ref=", "outfile=", "no_model_fit"])
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
        elif opt == "--no_model_fit":
            model_fit = False
        elif opt == "-b":
            num_bins = max(int(arg), 1)
        else:
            usage()
            sys.exit(1)

    if infile == '' or ref == '':
        print("Please specify the training reads and its reference genome!")
        usage()
        sys.exit(1)

    # READ PRE-PROCESS AND UNALIGNED READS ANALYSIS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")
    out_maf = outfile + ".maf"

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

    del dic_reads

    # if maf file provided
    if maf_file != '':
        if out_maf == maf_file:
            out_maf = outfile + "_processed.maf"

        call("grep '^s ' " + maf_file + " > " + out_maf, shell=True)

        # get best hit and unaligned reads
        unaligned_length = list(get_besthit.besthit_and_unaligned(in_fasta, out_maf, outfile))

    # if maf file not provided
    else:
        # Alignment
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
        call("lastdb ref_genome " + ref, shell=True)
        call("lastal -a 1 ref_genome " + in_fasta + " | grep '^s ' > " + out_maf, shell=True)

        # get best hit and unaligned reads
        unaligned_length = list(get_besthit.besthit_and_unaligned(in_fasta, out_maf, outfile))

    # ALIGNED READS ANALYSIS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
    num_aligned = align.head_align_tail(outfile, num_bins)

    # Length distribution of unaligned reads
    out1 = open(outfile + "_unaligned_length_ecdf", 'w')

    num_unaligned = len(unaligned_length)
    if num_unaligned != 0:
        max_length = max(unaligned_length)
        hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length + 50, 50),
                                                          density=True)
        cdf = numpy.cumsum(hist_unaligned * 50)
        out1.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
        out1.write("bin\t0-" + str(max_length) + '\n')
        for i in xrange(len(cdf)):
            out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    else:
        out1.write("Aligned / Unaligned ratio:\t100%\n")

    out1.close()
    del unaligned_length

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        path = sys.argv[0].split("/")
        r_path = '/'.join(path[:-1]) + '/' + "model_fitting.R"
        if os.path.isfile(r_path):
            call("R CMD BATCH '--args prefix=\"" + outfile + "\"' " + r_path, shell=True)
        else:
            sys.stderr.write("Could not find 'model_fitting.R' in ../src/\n" +
                  "Make sure you copied the whole source files from Github.")

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main(sys.argv[1:])

