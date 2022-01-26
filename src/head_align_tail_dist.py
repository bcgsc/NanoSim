#!/usr/bin/env python

from __future__ import with_statement
import numpy
import pysam
import sys
import joblib
from time import strftime
from sklearn.neighbors import KernelDensity
from file_handler import gzopen as open


def kde2d(x, y):
    x = numpy.array(x)
    y = numpy.array(y)
    xy = numpy.vstack([x, y])
    d = xy.shape[0]
    n = xy.shape[1]
    bw = (n * (d + 2) / 4.) ** (-1. / (d + 4))  # silverman
    kde_2d = KernelDensity(bandwidth=bw).fit(xy.T)

    return kde_2d


def edge_checker(rstart, rend, ref_length, ref_edge_max_dist=400, query_min_aln_len=100):
    is_edge = [False, False]
    if rend - rstart >= query_min_aln_len:
        if rend >= ref_length - 1 - ref_edge_max_dist:
            # read was aligned to reference end
            is_edge[1] = True
        elif rstart <= ref_edge_max_dist:
            # read was aligned to reference start
            is_edge[0] = True

    return is_edge


def get_head_tail(cigar_string, reverse_flag):
    head_info = cigar_string[0]
    tail_info = cigar_string[-1]

    if head_info[0] in (4, 5):  # soft clip is 4, hard clip is 5
        head = head_info[1]
    else:
        head = 0

    if tail_info[0] in (4, 5):
        tail = tail_info[1]
    else:
        tail = 0

    if reverse_flag:
        return tail, head
    else:
        return head, tail


def head_align_tail(prefix, alnm_ext, mode):
    alnm_file_prefix = prefix
    if mode == "transcriptome":
        total_ref_length = []  # total length of reference
        prefix = prefix[:-14]

        dict_genome_alnm_info = {}
        in_sam_file_genome = pysam.AlignmentFile(prefix + "_genome_primary.bam", 'rb')
        for alnm_temp in in_sam_file_genome.fetch(until_eof=True):
            read_temp = alnm_temp.query_name
            head_temp, tail_temp = get_head_tail(alnm_temp.cigartuples, alnm_temp.is_reverse)
            if read_temp not in dict_genome_alnm_info:
                dict_genome_alnm_info[read_temp] = (head_temp, tail_temp)
            else:
                h, t = dict_genome_alnm_info[read_temp]
                dict_genome_alnm_info[read_temp] = (min(head_temp, h), min(tail_temp, t))

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

                '''
                tail = int(query[5])-int(query[2])-int(query[3])
                ratio = float(query[3])/float(query[5])
                out1.write(query[5] + '\n')
                out2.write(query[3] + '\n')
                out3.write(query[2] + '\n')
                out4.write(str(aligned_ref) + '\n')
                out5.write(str(ht) + '\n')
                out6.write(str(ratio) + '\n')
                out7.write(str(tail) + '\n')
                '''
    else:
        last_read = ''
        last_ref = ''
        aligned_ref = 0
        in_sam_file = pysam.AlignmentFile(alnm_file_prefix + "_primary.bam", 'rb')

        # extract the lengths of all reference sequences
        dict_ref_len = {}
        for info in in_sam_file.header['SQ']:
            dict_ref_len[info['SN']] = info['LN']

        for alnm in in_sam_file.fetch(until_eof=True):
            read = alnm.query_name
            ref = alnm.reference_name
            if mode == "transcriptome":
                total_ref = dict_ref_len[ref]
                total_ref_length.append(total_ref)

            if read == last_read:
                # head and tail of a read with split alignments are considered together
                # aligned regions of a read with split alignments are considered separately
                # aligned regions of a circular read are considered are considered together
                if mode == "transcriptome" and (read in dict_genome_alnm_info):
                    head_g, tail_g = dict_genome_alnm_info[read]
                    head_t, tail_t = get_head_tail(alnm.cigartuples, alnm.is_reverse)
                    head_new = min(head_g, head_t)
                    tail_new = min(tail_g, tail_t)
                else:
                    head_new, tail_new = get_head_tail(alnm.cigartuples, alnm.is_reverse)
                head = min(head, head_new)
                tail = min(tail, tail_new)
                read_len_total = max(read_len_total, alnm.infer_read_length())
                if mode != "transcriptome":
                    is_edge = edge_checker(alnm.reference_start, alnm.reference_end, dict_ref_len[ref])

                # Check for "circular" reads
                if mode != "transcriptome" and ref == last_ref and \
                        ((last_is_edge[0] and is_edge[1]) or (last_is_edge[1] and is_edge[0])):
                    aligned_ref += alnm.reference_length
                    middle += alnm.query_alignment_length
                else:
                    aligned_ref_length.append(aligned_ref)
                    out2.write(str(last_read) + '\t' + str(middle) + '\n')
                    out4.write(str(last_read) + '\t' + str(aligned_ref) + '\n')
                    aligned_ref = alnm.reference_length
                    middle = alnm.query_alignment_length
                last_ref = ref

            else:
                if aligned_ref != 0:
                    aligned_ref_length.append(aligned_ref)
                    total_length.append(read_len_total)

                    # ratio aligned part over total length of the read
                    ratio = float(middle) / read_len_total
                    ht = head + tail
                    ht_length.append(ht)
                    if head != 0:
                        r = float(head) / ht
                        head_vs_ht_ratio.append(r)

                    out1.write(str(last_read) + '\t' + str(read_len_total) + '\n')
                    out2.write(str(last_read) + '\t' + str(middle) + '\n')
                    out3.write(str(last_read) + '\t' + str(head) + '\n')
                    out4.write(str(last_read) + '\t' + str(aligned_ref) + '\n')
                    out5.write(str(last_read) + '\t' + str(ht) + '\n')
                    out6.write(str(last_read) + '\t' + str(ratio) + '\n')
                    out7.write(str(last_read) + '\t' + str(tail) + '\n')

                    parsed += 1
                    if (parsed) % 1000 == 0:
                        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads processed >> " +
                                         str(parsed + 1) + "\r")
                        # +1 is just to ignore the zero index by python
                        sys.stdout.flush()

                last_read = read
                aligned_ref = alnm.reference_length
                read_len_total = alnm.infer_read_length()
                middle = alnm.query_alignment_length
                if mode == "transcriptome" and (read in dict_genome_alnm_info):
                    head_g = min(dict_genome_alnm_info[read][0])
                    tail_g = min(dict_genome_alnm_info[read][1])
                    head_t, tail_t = get_head_tail(alnm.cigartuples, alnm.is_reverse)
                    head = min(head_g, head_t)
                    tail = min(tail_g, tail_t)
                else:
                    head, tail = get_head_tail(alnm.cigartuples, alnm.is_reverse)
                if mode != "transcriptome":
                    last_is_edge = edge_checker(alnm.reference_start, alnm.reference_end, dict_ref_len[ref])
                last_ref = ref

        aligned_ref_length.append(aligned_ref)
        total_length.append(read_len_total)
        ratio = float(middle) / read_len_total
        ht = head + tail
        ht_length.append(ht)
        if ht != 0:
            r = float(head) / ht
            head_vs_ht_ratio.append(r)

        out1.write(str(last_read) + '\t' + str(read_len_total) + '\n')
        out2.write(str(last_read) + '\t' + str(middle) + '\n')
        out3.write(str(last_read) + '\t' + str(head) + '\n')
        out4.write(str(last_read) + '\t' + str(aligned_ref) + '\n')
        out5.write(str(last_read) + '\t' + str(ht) + '\n')
        out6.write(str(last_read) + '\t' + str(ratio) + '\n')
        out7.write(str(last_read) + '\t' + str(tail) + '\n')

    out1.close()
    out2.close()
    out3.close()
    out4.close()
    out5.close()
    out6.close()
    out7.close()

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing KDE\n")
    if mode == "transcriptome":
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing 2D KDE for transcriptome ref length\n")
        sys.stdout.flush()
        kde_2d = kde2d(total_ref_length, aligned_ref_length)
        joblib.dump(kde_2d, prefix + '_aligned_region_2d.pkl')

    aligned_length = numpy.array(aligned_ref_length)  # Aligned length of the reference, which is error-free
    total_length = numpy.array(total_length)
    ht_length = numpy.array(ht_length)
    head_vs_ht_ratio = numpy.array(head_vs_ht_ratio)

    # Kernel density function for the length of aligned regions
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing KDE for aligned region\n")
    sys.stdout.flush()
    aligned_2d = numpy.array(aligned_length)[:, numpy.newaxis]
    kde_aligned = KernelDensity(bandwidth=10).fit(aligned_2d)
    joblib.dump(kde_aligned, prefix + '_aligned_region.pkl')

    # Kernel density function for the length of aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing KDE for aligned reads\n")
    sys.stdout.flush()
    total_2d = numpy.array(total_length)[:, numpy.newaxis]
    kde_total = KernelDensity(bandwidth=10).fit(total_2d)
    joblib.dump(kde_total, prefix + '_aligned_reads.pkl')

    # Kernel density function for the length of ht
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing KDE for unaligned length\n")
    sys.stdout.flush()
    ht_log = numpy.log10(ht_length + 1)
    ht_log_2d = ht_log[:, numpy.newaxis]
    kde_ht = KernelDensity(bandwidth=0.01).fit(ht_log_2d)
    joblib.dump(kde_ht, prefix + '_ht_length.pkl')

    # Kernel density function for the head/total ratio
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Computing KDE for ht ratio\n")
    sys.stdout.flush()
    head_vs_ht_ratio_2d = head_vs_ht_ratio[:, numpy.newaxis]
    kde_ht_ratio = KernelDensity(bandwidth=0.01).fit(head_vs_ht_ratio_2d)
    joblib.dump(kde_ht_ratio, prefix + '_ht_ratio.pkl')

    num_aligned = len(total_length)
    return num_aligned
