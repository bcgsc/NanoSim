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
                    "-i : training ONT real reads, must be fasta files with extention '.fasta'\n" \
                    "-r : reference genome of the training reads\n" \
                    "-o : The prefix of output file, default = 'training'\n"

    sys.stderr.write(usage_message)


def main(argv):
    # Parse input and output files
    infile = ''
    outfile = 'training'
    ref = ''
    try:
        opts, args = getopt.getopt(argv, "hi:r:o:", ["infile=", "ref=", "outfile="])
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
        elif opt in ("-o", "--outfile"):
            outfile = arg
        else:
            usage()
            sys.exit(1)

    if infile == '' or ref == '':
        print("Please specify the training reads and its reference genome!")
        usage()
        sys.exit(1)

    # ALIGNMENT
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
    # call("lastdb ref_genome " + ref, shell=True)
    out_maf = outfile + ".maf"
    # call("lastal -a 1 ref_genome " + infile + " | grep '^s ' > " + out_maf, shell=True)
    '''
    # LENGTH DISTRIBUTION
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": length distribution analysis\n")
    # Get besthit of alignment
    unaligned_length = get_besthit.besthit_and_unaligned(infile, outfile)

    # Generate the histogram of aligned reads
    num_aligned = align.head_align_tail(outfile)

    # Length distribution of unaligned reads
    out1 = open(outfile + "_unaligned_length_ecdf", 'w')

    hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, 50001, 50), density=True)
    cdf = numpy.cumsum(hist_unaligned * 50)
    num_unaligned = len(unaligned_length)

    out1.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
    out1.write("bin\t0-50000\n")
    for i in xrange(len(cdf)):
        out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    out1.close()

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(outfile)
    '''
    path = sys.argv[0].split("/")
    R_path = '/'.join(path[:-1]) + '/' + "model_fitting.R"
    if os.path.isfile(R_path):
        call("Rscript " + R_path + " " + outfile, shell=True)
    else:
        print("Could not find 'model_fitting.R' in ../src/\nMake sure you copied the whole source files from Github.")
if __name__ == "__main__":
    main(sys.argv[1:])

