#!/usr/bin/env python
"""
Created on Sept 4th, 2015

@author: Chen Yang

This script transforms pairwise format alignment to maf format.

"""


import sys
import getopt
from file_handler import gzopen as open


def main(argv):
    # Parse input and output files
    infile = ''
    outfile = ''
    try:
        opts, args = getopt.getopt(argv, "hi:o:", ["ifile=", "ofile="])
    except getopt.GetoptError:
        print('pairwise2maf.py -i <inputfile> -o <outputfile>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print('test.py -i <inputfile>')
            sys.exit()
        elif opt in ("-i", "--ifile"):
            infile = arg
        elif opt in ("-o", "--ofile"):
            outfile = arg

    out = open(outfile, 'w')

    with open(infile, 'r') as f:
        for line in f:
            query = line.strip().split()
            qseq = next(f).strip()
            pair = next(f)
            rseq = next(f).strip()

            qid = query[0]
            if query[1] == "0":
                direction = "+"
            elif query[1] == "16":
                direction = "-"
            else:
                continue

            rid = query[2]
            rstart = int(query[3]) - 1
            cigar = query[5]

            parse_cigar = cigar.split("S")
            if len(parse_cigar) > 2:
                qstart = int(parse_cigar[0])
                qend = int(parse_cigar[1].split("M")[-1])
            elif len(parse_cigar) == 2:
                if parse_cigar[1] != "":
                    qstart = int(parse_cigar[0])
                    qend = 0
                else:
                    qstart = 0
                    qend = int(parse_cigar[0].split("M")[-1])
            else:
                qend = 0
                qstart = 0

            rseq = rseq[qstart: len(rseq) - qend]
            tmp = rseq.split("-")
            ralign = sum(len(x) for x in tmp)

            qseq = qseq[qstart: len(qseq) - qend]
            tmp = qseq.split("-")
            qalign = sum(len(x) for x in tmp)

            out.write("s " + rid + " " + str(rstart) + " " + str(ralign) + " + * " + rseq + '\n')
            out.write("s " + qid + " " + str(qstart) + " " + str(qalign) + " " + direction + " " +
                      str(qalign + qstart + qend) + " " + qseq + '\n')


if __name__ == "__main__":
    main(sys.argv[1:])
