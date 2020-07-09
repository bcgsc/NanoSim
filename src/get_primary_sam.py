#!/usr/bin/env python

from __future__ import with_statement
import pysam
import numpy
import re
import joblib
from sklearn.neighbors import KernelDensity


def cigar_parser(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    qstart = int(match[0][0]) if match[0][1] in ('S', 'H') else 0
    qlen = 0
    rlen = 0

    for item in match:
        if item[1] == 'M':
            qlen += int(item[0])
            rlen += int(item[0])
        elif item[1] == 'I':
            qlen += int(item[0])
        elif item[1] == 'D':
            rlen += int(item[0])
    qend = qstart + qlen
    return qstart, qend, qlen, rlen


def not_overlap(interval, interval_lst, interval_name=None, interval_name_list=None, overlap_base=10):
    # interval: (start, end)
    # overlap_base: a relaxed threshold reserved for a small stretch of bases shared by both alignments
    for i in range(len(interval_lst)):
        if interval[0] < interval_lst[i][1] - overlap_base and interval[1] - overlap_base > interval_lst[i][0]:
            if interval_name is None or (interval_name is not None and interval_name == interval_name_list[i]):
                return False
    return True


def primary_and_unaligned(sam_alnm_file, prefix):
    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=True)

    unaligned_len = []
    pos_strand = 0
    num_aligned = 0

    for aln in in_sam_file.fetch(until_eof=True):
        if not aln.is_unmapped and not aln.is_secondary and not aln.is_supplementary:
            num_aligned += 1
            out_sam_file.write(aln)
            if aln.flag == 0:
                pos_strand += 1
        elif aln.is_unmapped:
            unaligned_len.append(aln.query_length)

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    return unaligned_len, strandness


def primary_and_unaligned_circular(sam_alnm_file, prefix):
    """
    Function to extract alignments of reads at extremities of circular genome references
    :param sam_alnm_file: path of input SAM file
    :param prefix: prefix of output SAM file
    :outputs: a sam file that contains: 1) primary alignments for all reads that can be aligned and 2)
              split alignments for reads that are aligned to both ends of a circular genome,
              A pickle KDE file of the gap sizes between alignment fragments
    :returns: an unaligned_len list, and strandness information
    """

    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=True)
    gap_length = []
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0

    # extract the lengths of all reference sequences
    ref_lengths = dict()
    for info in in_sam_file.header['SQ']:
        ref_lengths[info['SN']] = info['LN']

    # parse alignment file, assuming that all alignments of a read are grouped together
    for aln in in_sam_file.fetch(until_eof=True):
        if aln.is_unmapped:
            unaligned_len.append(aln.query_length)

        elif not aln.is_secondary and not aln.is_supplementary:
            # a primary alignment
            num_aligned += 1

            # Define a list to store chimeric reads, each item is an interval (query_start, query_end)
            primary_direction = '+' if not aln.is_reverse else '-'
            NM_tag = int(aln.get_tag('NM'))
            primary_qstart = aln.query_alignment_start

            compatible_list = [{"query": [(aln.query_alignment_start, aln.query_alignment_end)],
                                "ref": [(aln.reference_start, aln.reference_end)],
                                "score": aln.query_alignment_length - NM_tag,
                                "rname": [aln.reference_name],
                                "direction": [primary_direction]}]

            supplementary_to_be_added = []
            try:
                supplementary_aln_list = aln.get_tag('SA').split(';')

                for supp_aln in supplementary_aln_list[:-1]:
                    ref_name, ref_start, direction, cigar, _ , NM_tag = supp_aln.split(',')
                    ref_start = int(ref_start) - 1
                    NM_tag = int(NM_tag)
                    qstart, qend, qlen, rlen = cigar_parser(cigar)

                    added = False
                    for seg in compatible_list:
                        if not_overlap((qstart, qend), seg["query"]) and \
                                not_overlap((ref_start, ref_start + rlen), seg["ref"], ref_name, seg['rname']):
                            seg["query"].append((qstart, qend))
                            seg["ref"].append((ref_start, ref_start + rlen))
                            seg["score"] += (qlen - NM_tag)
                            seg["rname"].append(ref_name)
                            seg["direction"].append(direction)
                            added = True

                    if not added:
                        compatible_list.append({"query": [(qstart, qend)],
                                                "ref": [(ref_start, ref_start + rlen)],
                                                "score": qlen - NM_tag,
                                                "rname": [ref_name],
                                                "direction": [direction]})

                max_score = max([x["score"] for x in compatible_list])
                for seg in compatible_list:
                    if seg["score"] == max_score:
                        idx = [i[0] for i in sorted(enumerate(seg["query"]), key=lambda x:x[1])]
                        seg["query"].sort()
                        seg["ref"] = [seg["ref"][x] for x in idx]
                        seg["rname"] = [seg["rname"][x] for x in idx]

                        dir_added = False
                        for i in range(len(seg["query"])):
                            interval = seg["query"][i]
                            if i > 0:
                                gap = max(0, interval[0] - seg["query"][i - 1][1])  # Change negative gaps size to 0
                                gap_length.append(gap)
                            if interval[0] == primary_qstart:
                                dir_added = True
                                if primary_direction == '+':
                                    pos_strand += 1
                                out_sam_file.write(aln)
                            else:
                                supplementary_to_be_added.append((seg['rname'][i], seg['query'][i][0],
                                                                 seg['ref'][i][0]))
                        if not dir_added:
                            if seg["direction"][0] == '+':
                                pos_strand += 1
                        break

            except KeyError:
                out_sam_file.write(aln)
                if primary_direction == '+':
                    pos_strand += 1

        else:
            qstart, _, _, _ = cigar_parser(aln.cigarstring)
            if (aln.reference_name, qstart, aln.reference_start) in supplementary_to_be_added:
                out_sam_file.write(aln)

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    # Compute the KDE of gaps
    gap_length = numpy.array(gap_length)
    gap_log = numpy.log10(gap_length + 1)
    gap_log_2d = gap_log[:, numpy.newaxis]
    kde_gap = KernelDensity(bandwidth=0.01).fit(gap_log_2d)
    joblib.dump(kde_gap, prefix + '_gap_length.pkl')

    return unaligned_len, strandness

