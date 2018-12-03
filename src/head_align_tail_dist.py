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
import HTSeq
from sklearn.neighbors import KernelDensity
from math import log, ceil
from sklearn.externals import joblib

try:
    from six.moves import xrange
except ImportError:
    pass


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


def get_head_tail(cigar_string):

    head_info = cigar_string[0]
    tail_info = cigar_string[-1]

    if head_info.type == "S" or head_info.type == "H":
        head = head_info.size
    else:
        head = 0

    if tail_info.type == "S" or tail_info.type == "H":
        tail = tail_info.size
    else:
        tail = 0

    return head, tail


def head_align_tail(prefix, alnm_ftype):
    out1 = open(prefix + "_total.txt", 'w')
    out2 = open(prefix + "_middle.txt", 'w')
    out3 = open(prefix + "_head.txt", 'w')
    out4 = open(prefix + "_middle_ref.txt", 'w')
    out5 = open(prefix + "_ht.txt", 'w')
    out6 = open(prefix + "_ratio.txt", 'w')
    out7 = open(prefix + "_tail.txt", 'w')

    aligned_length = []
    total_length = []
    ht_length = []
    head_vs_ht_ratio = []

    if alnm_ftype == "maf":
        besthit_out = prefix + "_besthit.maf"
        with open(besthit_out, 'r') as f:
            for line in f:
                ref = line.strip().split()
                aligned_ref = int(ref[3])
                aligned_length.append(aligned_ref)
                query = next(f).strip().split()
                head = int(query[2])
                middle = int(query[3])
                tail = int(query[5])-int(query[2])-int(query[3])
                total_length.append(int(query[5]))
                ht = int(query[5])-int(query[3])
                ht_length.append(ht)
                ratio = float(query[3])/float(query[5])

                if ht != 0:
                    r = float(head) / ht
                    head_vs_ht_ratio.append(r)
                out1.write(query[5] + '\n')
                out2.write(query[3] + '\n')
                out3.write(query[2] + '\n')
                out4.write(str(aligned_ref) + '\n')
                out5.write(str(ht) + '\n')
                out6.write(str(ratio) + '\n')
                out7.write(str(tail) + '\n')
    else:
        sam_reader = HTSeq.SAM_Reader
        alnm_file_sam = prefix + "_primary.sam"
        alignments = sam_reader(alnm_file_sam)
        for alnm in alignments:
            ref = alnm.iv.chrom
            aligned_ref = alnm.iv.length
            aligned_length.append(aligned_ref)

            read_len_total = len(alnm.read.seq)
            total_length.append(read_len_total)
            head, tail = get_head_tail(alnm.cigar)
            middle = read_len_total - head - tail

            # ratio aligned part over total length of the read
            ratio = float(middle) / read_len_total
            ht = head + tail
            ht_length.append(ht)
            if head != 0:
                r = float(head) / ht
                head_vs_ht_ratio.append(r)
            out1.write(str(read_len_total) + '\n')
            out2.write(str(middle) + '\n')
            out3.write(str(head) + '\n')
            out4.write(str(aligned_ref) + '\n')
            out5.write(str(ht) + '\n')
            out6.write(str(ratio) + '\n')
            out7.write(str(tail) + '\n')

    out1.close()
    out2.close()
    out3.close()
    out4.close()
    out5.close()
    out6.close()
    out7.close()

    aligned_length = numpy.array(aligned_length)
    total_length = numpy.array(total_length)
    ht_length = numpy.array(ht_length)
    head_vs_ht_ratio = numpy.array(head_vs_ht_ratio)

    # Kernel density function for the length of aligned regions
    aligned_2d = numpy.array(aligned_length)[:, numpy.newaxis]
    kde_aligned = KernelDensity(bandwidth=10).fit(aligned_2d)
    joblib.dump(kde_aligned, prefix + '_aligned_region.pkl')

    # Kernel density function for the length of aligned reads
    total_2d = numpy.array(total_length)[:, numpy.newaxis]
    kde_total = KernelDensity(bandwidth=10).fit(total_2d)
    joblib.dump(kde_total, prefix + '_aligned_reads.pkl')

    # Kernel density function for the length of ht
    ht_log = numpy.log10(ht_length + 1)
    ht_log_2d = ht_log[:, numpy.newaxis]
    kde_ht = KernelDensity(bandwidth=0.01).fit(ht_log_2d)
    joblib.dump(kde_ht, prefix + '_ht_length.pkl')

    # Kernel density function for the head/total ratio
    head_vs_ht_ratio_2d = head_vs_ht_ratio[:, numpy.newaxis]
    kde_ht_ratio = KernelDensity(bandwidth=0.01).fit(head_vs_ht_ratio_2d)
    joblib.dump(kde_ht_ratio, prefix + '_ht_ratio.pkl')

    num_aligned = len(total_length)
    return num_aligned
