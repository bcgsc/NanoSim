#!/usr/bin/env python
"""
Written by Chen Yang on Mar 25th, 2015
To get the length of head, aligned, and tail regions of an alignment.

Major change in Apr 22nd

Updated in Nov 25th
"""

from __future__ import with_statement
import sys
import getopt
import numpy


def head_align_tail(outfile):
    out1 = open(outfile + '_aligned_length_ecdf', 'w')
    out2 = open(outfile + '_aligned_reads_ecdf', 'w')
    out3 = open(outfile + '_ht_ratio', 'w')
    out4 = open(outfile + "_align_ratio", 'w')
    aligned = []
    total = []
    ht_ratio = {
        (0, 10): [], (10, 20): [], (20, 30): [], (30, 40): [], (40, 50): [], (50, 100): [], (100, 300): [], (300, 1000): [],
        (1000, 2000): [], (2000, 3000): [], (3000, 5000): [], (5000, 7500): [], (7500, 10000): [], (10000, 50000): []}
    align_ratio = {
        (0, 1000): [], (1000, 2000): [], (2000, 3000): [], (3000, 4000): [], (4000, 5000): [], (5000, 6000): [],
        (6000, 7000): [], (7000, 8000): [], (8000, 9000): [], (9000, 10000): [], (10000, 11000): [], (11000, 12000): [],
        (12000, 13000): [], (13000, 15000): [], (15000, 20000): [], (20000, 25000): [], (25000, 50000): []}

    besthit_out = outfile + "_besthit.maf"
    with open(besthit_out, 'r') as f:
        for line in f:
            new = line.strip().split()
            aligned_ref = int(new[3])
            aligned.append(aligned_ref)
            new = next(f).strip().split()
            head = int(new[2])
            # tail = int(new[5])-int(new[2])-int(new[3])
            total.append(int(new[5]))
            ht = int(new[5])-int(new[3])
            ratio = float(new[3])/float(new[5])
            for key in align_ratio.keys():
                if key[0] <= int(new[3]) < key[1]:
                    align_ratio[key].append(ratio)
            for k in ht_ratio.keys():
                if k[0] <= ht < k[1]:
                    if ht != 0:
                        r = float(head) / ht
                        ht_ratio[k].append(r)
                    break

    # ecdf of align ratio
    align_cum = {(0, 1000): [], (1000, 2000): [], (2000, 3000): [], (3000, 4000): [], (4000, 5000): [],
                 (5000, 6000): [], (6000, 7000): [], (7000, 8000): [], (8000, 9000): [], (9000, 10000): [],
                 (10000, 11000): [], (11000, 12000): [], (12000, 13000): [], (13000, 15000): [], (15000, 20000): [],
                 (20000, 25000): [], (25000, 50000): []}
    for key, value in align_ratio.items():
        if value:
            hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
            cdf = numpy.cumsum(hist_ratio * 0.001)
        else:
            cdf = [0] * 1000
        align_cum[key] = cdf

    out4.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(align_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out4.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t")
        for key in sorted(align_cum.keys()):
            out4.write(str(align_cum[key][i]) + "\t")
        out4.write("\n")

    # ecdf of head/total ratio
    ht_cum = {(0, 10): [], (10, 20): [], (20, 30): [], (30, 40): [], (40, 50): [], (50, 100): [], (100, 300): [],
              (300, 1000): [], (1000, 2000): [], (2000, 3000): [], (3000, 5000): [], (5000, 7500): [], (7500, 10000): [],
              (10000, 50000): []}
    for key, value in ht_ratio.items():
        if value:
            hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
            cdf = numpy.cumsum(hist_ratio * 0.001)
        else:
            cdf = [0] * 1000
        ht_cum[key] = cdf

    out3.write("bins\t" + "\t".join("%s-%s" % tup for tup in sorted(ht_cum.keys())) + "\n")
    for i in xrange(len(cdf)):
        out3.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t")
        for key in sorted(ht_cum.keys()):
            out3.write(str(ht_cum[key][i]) + "\t")
        out3.write("\n")

    # ecdf of length of aligned regions
    hist_aligned, bin_edges = numpy.histogram(aligned, bins=numpy.arange(0, 50001, 50), density=True)
    cdf = numpy.cumsum(hist_aligned * 50)
    out1.write("bin\t0-50000\n")
    for i in xrange(len(cdf)):
        out1.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t" + str(cdf[i]) + '\n')
    num_aligned = len(aligned)

    # ecdf of length of aligned reads
    hist_reads, bin_edges = numpy.histogram(total, bins=numpy.arange(0, 50001, 50), density=True)
    cdf = numpy.cumsum(hist_reads * 50)
    out2.write("bin\t0-50000\n")
    for i in xrange(len(cdf)):
        out2.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t" + str(cdf[i]) + '\n')

    out1.close()
    out2.close()
    out3.close()
    out4.close()
    return num_aligned
