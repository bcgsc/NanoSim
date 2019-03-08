#!/usr/bin/env python
"""
Written by Chen Yang on Mar 25th, 2015
To get the length of head, aligned, and tail regions of an alignment.

Major change in Apr 22nd

Updated in Nov 25th
"""

from __future__ import with_statement
import numpy
import HTSeq
from sklearn.neighbors import KernelDensity
from math import log, ceil
from sklearn.externals import joblib

try:
    from six.moves import xrange
except ImportError:
    pass


def kde2d(x, y):
    x = numpy.array(x)
    y = numpy.array(y)
    xy = numpy.vstack([x, y])
    d = xy.shape[0]
    n = xy.shape[1]
    bw = (n * (d + 2) / 4.) ** (-1. / (d + 4))  # silverman
    kde_2d = KernelDensity(bandwidth=bw).fit(xy.T)
    # xmin = x.min()
    # xmax = x.max()
    # ymin = y.min()
    # ymax = y.max()
    # X, Y = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    # positions = np.vstack([X.ravel(), Y.ravel()])
    # Z = np.reshape(np.exp(kde_2d.score_samples(positions.T)), X.shape)

    return kde_2d


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


def head_align_tail(*args):
    '''
    out1 = open(prefix + "_total.txt", 'w')
    out2 = open(prefix + "_middle.txt", 'w')
    out3 = open(prefix + "_head.txt", 'w')
    out4 = open(prefix + "_middle_ref.txt", 'w')
    out5 = open(prefix + "_ht.txt", 'w')
    out6 = open(prefix + "_ratio.txt", 'w')
    out7 = open(prefix + "_tail.txt", 'w')
    '''
    prefix = args[0]
    alnm_ftype = args[1]
    mode = args[2]
    if len(args) == 4:
        dict_ref_len = args[3]


    if mode == "transcriptome":
        x = []  # total length of reference
        y = []  # aligned length on reference

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
                if mode == "transcriptome":
                    total_ref = int(ref[5])
                    x.append(total_ref)
                    y.append(aligned_ref)
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
                '''
                out1.write(query[5] + '\n')
                out2.write(query[3] + '\n')
                out3.write(query[2] + '\n')
                out4.write(str(aligned_ref) + '\n')
                out5.write(str(ht) + '\n')
                out6.write(str(ratio) + '\n')
                out7.write(str(tail) + '\n')
                '''
    else:
        sam_reader = HTSeq.SAM_Reader
        alnm_file_sam = prefix + "_primary.sam"
        alignments = sam_reader(alnm_file_sam)
        for alnm in alignments:
            ref = alnm.iv.chrom
            aligned_ref = alnm.iv.length
            aligned_length.append(aligned_ref)
            if mode == "transcriptome":
                total_ref = dict_ref_len[ref]
                x.append(total_ref)
                y.append(aligned_ref)

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
            '''
            out1.write(str(read_len_total) + '\n')
            out2.write(str(middle) + '\n')
            out3.write(str(alnm.read.name) + '\t' + str(head) + '\n')
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
    '''
    if mode == "transcriptome":
        kde_2d = kde2d(x, y)
        joblib.dump(kde_2d, prefix + '_aligned_region_2d.pkl')

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
