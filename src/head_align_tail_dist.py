#!/usr/bin/env python

from __future__ import with_statement
import numpy
import HTSeq
import sys
import joblib
from time import strftime
from sklearn.neighbors import KernelDensity


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
    prefix = args[0]
    alnm_ext = args[1]
    mode = args[2]
    alnm_file_prefix = prefix
    if mode == "transcriptome":
        total_ref_length = []  # total length of reference
        dict_ref_len = args[3]
        prefix = prefix[:-14]

    aligned_ref_length = []  # aligned length of reference
    total_length = []
    ht_length = []
    head_vs_ht_ratio = []

    # For debugging
    out1 = open(prefix + "_total.txt", 'w')
    out2 = open(prefix + "_middle.txt", 'w')
    out3 = open(prefix + "_head.txt", 'w')
    out4 = open(prefix + "_middle_ref.txt", 'w')
    out5 = open(prefix + "_ht.txt", 'w')
    out6 = open(prefix + "_ratio.txt", 'w')
    out7 = open(prefix + "_tail.txt", 'w')

    parsed = 0

    if alnm_ext == "maf":
        # TODO: circular reads
        alnm_file = alnm_file_prefix + "_besthit.maf"
        with open(alnm_file, 'r') as f:
            for line in f:
                ref = line.strip().split()
                aligned_ref = int(ref[3])
                if mode == "transcriptome":
                    total_ref = int(ref[5])  # Need further investigation
                    total_ref_length.append(total_ref)
                aligned_ref_length.append(aligned_ref)
                query = next(f).strip().split()
                head = int(query[2])
                total_length.append(int(query[5]))
                ht = int(query[5])-int(query[3])
                ht_length.append(ht)

                if ht != 0:
                    r = float(head) / ht
                    head_vs_ht_ratio.append(r)

                tail = int(query[5])-int(query[2])-int(query[3])
                ratio = float(query[3])/float(query[5])
                out1.write(query[5] + '\n')
                out2.write(query[3] + '\n')
                out3.write(query[2] + '\n')
                out4.write(str(aligned_ref) + '\n')
                out5.write(str(ht) + '\n')
                out6.write(str(ratio) + '\n')
                out7.write(str(tail) + '\n')
    else:
        sam_reader = HTSeq.SAM_Reader
        alnm_file = alnm_file_prefix + "_primary.sam"
        alignments = sam_reader(alnm_file)
        last_read = ''
        aligned_ref = 0
        for alnm in alignments:
            read = alnm.read.name
            ref = alnm.iv.chrom
            if mode == "transcriptome":
                total_ref = dict_ref_len[ref]
                total_ref_length.append(total_ref)
                aligned_ref_length.append(alnm.iv.length)
                read_len_total = len(alnm.read.seq)
                total_length.append(read_len_total)
                head, tail = get_head_tail(alnm.cigar)

            else:
                if read == last_read:
                    flag = False
                    aligned_ref += alnm.iv.length
                    head_new, tail_new = get_head_tail(alnm.cigar)
                    head = min(head, head_new)
                    tail = min(tail, tail_new)
                else:
                    flag = True
                    if aligned_ref != 0:
                        aligned_ref_length.append(aligned_ref)
                        total_length.append(read_len_total)
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
                        out3.write(str(last_read) + '\t' + str(head) + '\n')
                        out4.write(str(aligned_ref) + '\n')
                        out5.write(str(ht) + '\n')
                        out6.write(str(ratio) + '\n')
                        out7.write(str(tail) + '\n')
                        parsed += 1
                        if (parsed) % 1000 == 0:
                            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads processed >> " +
                                             str(parsed + 1) + "\r")
                            # +1 is just to ignore the zero index by python
                            sys.stdout.flush()

                    last_read = alnm.read.name
                    aligned_ref = alnm.iv.length
                    read_len_total = len(alnm.read.seq)
                    head, tail = get_head_tail(alnm.cigar)

    if not flag:
        aligned_ref_length.append(aligned_ref)
        total_length.append(read_len_total)
        middle = read_len_total - head - tail
        ratio = float(middle) / read_len_total
        ht = head + tail
        ht_length.append(ht)
        if head != 0:
            r = float(head) / ht
            head_vs_ht_ratio.append(r)

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

    if mode == "transcriptome":
        kde_2d = kde2d(total_ref_length, aligned_ref_length)
        joblib.dump(kde_2d, prefix + '_aligned_region_2d.pkl')

    aligned_length = numpy.array(aligned_ref_length)  # Aligned length of the reference, which is error-free
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
