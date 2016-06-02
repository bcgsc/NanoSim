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


def flex_bins(num_of_bins, ratio_dict, num_of_reads):
    count_reads = num_of_reads / num_of_bins
    k_of_bin = 0
    k_of_ratio = 0
    ratio_keys = sorted(ratio_dict.keys())
    num_of_keys = len(ratio_keys)

    ratio_bins = {}
    while k_of_bin < num_of_bins:
        if k_of_ratio >= num_of_keys:
            break

        start = k_of_ratio
        count = len(ratio_dict[ratio_keys[k_of_ratio]])
        k_of_ratio += 1

        while k_of_ratio < num_of_keys:
            tmp_count = count + len(ratio_dict[ratio_keys[k_of_ratio]])
            if abs(tmp_count - count_reads) >= abs(count - count_reads):
                break
            else:
                count = tmp_count
                k_of_ratio += 1

        k = (ratio_keys[start] if start else 0,
             ratio_keys[k_of_ratio] if k_of_ratio < num_of_keys else ratio_keys[k_of_ratio - 1] + 1)
        ratio_bins[k] = []
        for i in xrange(start, k_of_ratio):
            ratio_bins[k].extend(ratio_dict[ratio_keys[i]])

        k_of_bin += 1

    if k_of_ratio < num_of_keys - 1:
        k = (ratio_keys[k_of_ratio], ratio_keys[num_of_keys - 1] + 1)
        ratio_bins[k] = []
        for i in xrange(k_of_ratio, num_of_keys - 1):
            ratio_bins[k].extend(ratio_dict[ratio_keys[i]])

    return ratio_bins


def head_align_tail(outfile, num_of_bins):
    out1 = open(outfile + '_aligned_length_ecdf', 'w')
    out2 = open(outfile + '_aligned_reads_ecdf', 'w')
    out3 = open(outfile + '_ht_ratio', 'w')
    out4 = open(outfile + "_align_ratio", 'w')

    '''
    out5 = open(outfile + "_total.txt", 'w')
    out6 = open(outfile + "_middle.txt", 'w')
    out7 = open(outfile + "_head.txt", 'w')
    out8 = open(outfile + "_middle_ref.txt", 'w')
    out9 = open(outfile + "_ht.txt", 'w')
    out10 = open(outfile + "_ratio.txt", 'w')
    out11 = open(outfile + "_tail.txt", 'w')
    '''

    aligned = []
    total = []
    ht_ratio = {}
    align_ratio = {}

    besthit_out = outfile + "_besthit.maf"
    with open(besthit_out, 'r') as f:
        for line in f:
            ref = line.strip().split()
            aligned_ref = int(ref[3])
            aligned.append(aligned_ref)
            query = next(f).strip().split()
            head = int(query[2])
            middle = int(query[3])
            tail = int(query[5])-int(query[2])-int(query[3])
            total.append(int(query[5]))
            ht = int(query[5])-int(query[3])
            ratio = float(query[3])/float(query[5])
            if middle in align_ratio:
                align_ratio[middle].append(ratio)
            else:
                align_ratio[middle] = [ratio]
            if ht != 0:
                r = float(head) / ht
                if ht in ht_ratio:
                    ht_ratio[ht].append(r)
                else:
                    ht_ratio[ht] = [r]
    '''
            out5.write(query[5] + '\n')
            out6.write(query[3] + '\n')
            out7.write(query[2] + '\n')
            out8.write(str(aligned_ref) + '\n')
            out9.write(str(ht) + '\n')
            out10.write(str(ratio) + '\n')
            out11.write(str(tail) + '\n')
    
    out5.close()
    out6.close()
    out7.close()
    out8.close()
    out9.close()
    out10.close()
    out11.close()
    '''

    max_length = max(total)

    # ecdf of length of aligned regions
    hist_aligned, bin_edges = numpy.histogram(aligned, bins=numpy.arange(0, max_length + 1, 1), density=True)
    cdf = numpy.cumsum(hist_aligned * 1)
    out1.write("bin\t0-" + str(max_length) + '\n')
    for i in xrange(len(cdf)):
        out1.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t" + str(cdf[i]) + '\n')
    num_aligned = len(aligned)

    # ecdf of length of aligned reads
    hist_reads, bin_edges = numpy.histogram(total, bins=numpy.arange(0, max_length + 1, 1), density=True)
    cdf = numpy.cumsum(hist_reads * 1)
    out2.write("bin\t0-" + str(max_length) + '\n')
    for i in xrange(len(cdf)):
        out2.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t" + str(cdf[i]) + '\n')

    # ecdf of head/total ratio
    # there needs to be at least one bin

    ht_ratio_bins = flex_bins(num_of_bins, ht_ratio, len(total))

    ht_cum = dict.fromkeys(ht_ratio_bins.keys(), [])
    for key, value in ht_ratio_bins.items():
        hist_ht, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ht * 0.001)
        ht_cum[key] = cdf

    out3.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(ht_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out3.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t")
        for key in sorted(ht_cum.keys()):
            out3.write(str(ht_cum[key][i]) + "\t")
        out3.write("\n")

    # ecdf of align ratio
    align_ratio_bins = flex_bins(num_of_bins, align_ratio, len(total))

    align_cum = dict.fromkeys(align_ratio_bins.keys(), [])
    for key, value in align_ratio_bins.items():
        hist_ratio, bin_edges = numpy.histogram(value, bins=numpy.arange(0, 1.001, 0.001), density=True)
        cdf = numpy.cumsum(hist_ratio * 0.001)
        align_cum[key] = cdf

    out4.write("bins\t" + '\t'.join("%s-%s" % tup for tup in sorted(align_cum.keys())) + '\n')
    for i in xrange(len(cdf)):
        out4.write(str(bin_edges[i]) + '-' + str(bin_edges[i+1]) + "\t")
        for key in sorted(align_cum.keys()):
            out4.write(str(align_cum[key][i]) + "\t")
        out4.write("\n")

    out1.close()
    out2.close()
    out3.close()
    out4.close()
    return num_aligned
