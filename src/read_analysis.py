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
import get_besthit_maf
import get_primary_sam
import besthit_to_histogram as error_model


# Usage information
def usage():
    usage_message = "./read_analysis.py <options>\n" \
                    "<options>: \n" \
                    "-h : print usage message\n" \
                    "-i : training ONT real reads, must be fasta files\n" \
                    "-r : reference genome of the training reads\n" \
                    "-a : Aligner to be used: minimap2 or LAST, default = 'minimap2\n" \
                    "-m : User can provide their own alignment file, with maf or sam extension, can be omitted\n" \
                    "-b : number of bins (for development), default = 20\n" \
                    "-t : number of threads for LAST alignment, default = 1\n" \
                    "-o : The prefix of output file, default = 'training'\n" \
                    "--no_model_fit : Skip the model fitting step\n"

    sys.stderr.write(usage_message)


def main(argv):
    # Parse input and output files
    infile = ''
    prefix = 'training'
    ref = ''
    aligner = ''
    alnm_file = ''
    model_fit = True
    num_threads = '1'
    num_bins = 20
    try:
        opts, args = getopt.getopt(argv, "hi:r:a:o:m:b:t:", ["infile=", "ref=", "prefix=", "no_model_fit"])
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
        elif opt == "-a":
            aligner = arg
        elif opt == "-m":
            alnm_file = arg
        elif opt in ("-o", "--prefix"):
            prefix = arg
        elif opt == "--no_model_fit":
            model_fit = False
        elif opt == "-b":
            num_bins = max(int(arg), 1)
        elif opt == "-t":
            num_threads = arg
        else:
            usage()
            sys.exit(1)

    if infile == '' or ref == '':
        print("Please specify the training reads and its reference genome!")
        usage()
        sys.exit(1)

    if aligner != '' and alnm_file != '':
        print("Please specify either an alignment file (-m ) OR an aligner to use for alignment (-a )")
        usage()
        sys.exit(1)

    # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")

    # Read pre-process
    in_fasta = prefix + "_processed.fasta"  # use the prefix of input fasta file for processed fasta file
    processed_fasta = open(in_fasta, 'w')
    dic_reads = {}
    with open(infile, 'r') as f:
        for line in f:
            if line[0] == '>':
                name = '-'.join(line.strip()[1:].split())
                dic_reads[name] = ""
            else:
                dic_reads[name] += line.strip()
    for k, v in dic_reads.items():
        processed_fasta.write('>' + k + '\n' + v + '\n')
    processed_fasta.close()

    del dic_reads

    # if an alignment file is provided
    if alnm_file != '':
        pre, file_ext = os.path.splitext(alnm_file)
        file_extension = file_ext[1:]
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing alignment file: " + file_extension + "\n")
        if file_extension == "maf":
            processed_maf = prefix + "_processed.maf"

            call("grep '^s ' " + alnm_file + " > " + processed_maf, shell=True)

            # get best hit and unaligned reads
            unaligned_length = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf, prefix)

        elif file_extension == "sam":
            #get the primary alignments and define unaligned reads.
            unaligned_length = get_primary_sam.primary_and_unaligned(alnm_file, prefix)
        else:
            print("Please specify an acceptable alignment format! (maf or sam)\n")
            usage()
            sys.exit(1)

    # if alignment file is not provided
    else:
        if aligner == "minimap2" or aligner == "": # Align with minimap2 by default
            file_extension = "sam"
            out_sam = prefix + ".sam"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2\n")
            call("minimap2 --cs -ax map-ont " + ref + " " + in_fasta + " > " + out_sam, shell=True)
            # get primary alignments and unaligned reads
            unaligned_length = get_primary_sam.primary_and_unaligned(out_sam, prefix)
        elif aligner == "LAST":
            file_extension = "maf"
            out_maf = prefix + ".maf"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
            call("lastdb ref_genome " + ref, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " | grep '^s ' > " + out_maf, shell=True)
            unaligned_length = get_besthit_maf.besthit_and_unaligned(in_fasta, out_maf, prefix)
        else:
            print("Please specify an acceptable aligner (minimap2 or LAST)\n")
            usage()
            sys.exit(1)


    # Aligned reads analysis
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
    num_aligned = align.head_align_tail(prefix, num_bins, file_extension)

    # Length distribution of unaligned reads
    out_unaligned_ecdf = open(prefix + "_unaligned_length_ecdf", 'w')

    num_unaligned = len(unaligned_length)
    if num_unaligned != 0:
        max_length = max(unaligned_length)
        hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length + 50, 50),
                                                          density=True)
        cdf = numpy.cumsum(hist_unaligned * 50)
        out_unaligned_ecdf.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
        out_unaligned_ecdf.write("bin\t0-" + str(max_length) + '\n')
        for i in xrange(len(cdf)):
            out_unaligned_ecdf.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
    else:
        out_unaligned_ecdf.write("Aligned / Unaligned ratio:\t100%\n")

    out_unaligned_ecdf.close()
    del unaligned_length

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    error_model.hist(prefix, file_extension)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        path = sys.argv[0].split("/")
        r_path = '/'.join(path[:-1]) + '/' + "model_fitting.R"
        if os.path.isfile(r_path):
            call("R CMD BATCH '--args prefix=\"" + prefix + "\"' " + r_path, shell=True)
        else:
            sys.stderr.write("Could not find 'model_fitting.R' in ../src/\n" +
                  "Make sure you copied the whole source files from Github.")

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main(sys.argv[1:])

