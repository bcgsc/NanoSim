#!/usr/bin/env python

from __future__ import with_statement
import pysam
import numpy
import re


def cigar_parser(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    qstart = int(match[0][0])
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


def not_overlap(interval, interval_lst, overlap_base=10):
    # interval: (start, end)
    # overlap_base: a relaxed threshold reserved for a small stretch of bases shared by both alignments
    for i in interval_lst:
        if interval[0] < i[1] - overlap_base and interval[1] - overlap_base > i[0]:
            return False
    return True


def edge_checker(rstart, rend, rname, ref_lengths, ref_edge_max_dist=400, query_min_aln_len=100):
    is_edge = [False, False]
    if rend - rstart >= query_min_aln_len:
        if rend >= ref_lengths[rname] - 1 - ref_edge_max_dist:
            # read was aligned to reference end
            is_edge[1] = True
        elif rstart <= ref_edge_max_dist:
            # read was aligned to reference start
            is_edge[0] = True

    return is_edge


def primary_and_unaligned(sam_alnm_file, prefix):
    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=False)

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


def primary_and_unaligned_circular(sam_alnm_file, prefix, ref_edge_max_dist=400, query_min_aln_len=100,
                                   overlap_base=10):
    """
    Function to extract alignments of reads at extremities of circular genome references
    :param sam_alnm_file: path of input SAM file
    :param prefix: prefix of output SAM file
    :param meta_list: a dictionary of metagenome information, including circularity
    :param ref_edge_max_dist: max distance of alignments from extremities of references
    :param query_min_aln_len: min length of split alignments
    :return: Outputs a sam file that contains: 1) primary alignments for all reads that can be aligned and 2)
             split alignments for reads that are aligned to both ends of a circular genome
             Return unaligned_len list, and strandness information
    """

    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=False)
    tmp = open("chimeric", 'w')
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
            continue

        if not aln.is_secondary and not aln.is_supplementary:
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
                                not_overlap((ref_start, ref_start + rlen), seg["ref"]):
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
                        out_info = "Query " + ';'.join(str(x[0]) + '-' + str(x[1]) for x in seg["query"])
                        out_info += '\t' + "Ref " + ';'.join(str(x[0]) + '-' + str(x[1]) for x in seg["ref"])
                        out_info += '\t' + "Direction " + ';'.join(x for x in seg["direction"]) + '\n'

                        tmp.write(aln.query_name + '\t' + out_info)

                '''
                edge_added = [False, False]
                added_strand = ''
                for seg in compatible_list:
                    if seg["score"] == max_score:
                        for i in range(len(seg['query'])):
                            edge_info = edge_checker(seg['ref'][i][0], seg['ref'][i][1], seg['rname'][i], ref_lengths)
                            if seg['query'][i][0] == primary_qstart:
                                out_sam_file.write(aln)
                                edge_added = edge_info
                                added_strand = seg['direction'][i]
                                if seg['direction'][i] == '+':
                                    pos_strand += 1
                                if not any(edge_added):
                                    break
                            elif (edge_info[0] and not edge_added[0]) or (edge_info[1] and not edge_added[1]):
                                if added_strand == '' or added_strand == seg['direction'][i]:
                                    if added_strand == '' and seg['direction'][i] == '+':
                                        pos_strand += 1
                                    added_strand = seg['direction'][i]
                                supplementary_to_be_added.append(seg['ref'][i][0])
                            else:
                                continue
                    else:
                        continue
                '''

            except KeyError:
                out_sam_file.write(aln)
                if primary_direction == '+':
                    pos_strand += 1

        else:
            continue
            '''
            if aln.reference_start in supplementary_to_be_added:
                out_sam_file.write(aln)
            '''

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    tmp.close()

    return unaligned_len, strandness

