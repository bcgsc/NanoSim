#!/usr/bin/env python

from __future__ import with_statement
import HTSeq
import pysam
import numpy


def primary_and_unaligned(sam_alnm_file, prefix):
    out_primary = open(prefix + "_primary.sam", 'w')
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
    sam_reader = HTSeq.SAM_Reader
    alignments = sam_reader(sam_alnm_file)
    for aln in alignments:
        if aln.aligned and (not aln.not_primary_alignment and not aln.supplementary):
            out_primary.write(aln.original_sam_line)
            num_aligned += 1
            if aln.flag == 0:
                pos_strand += 1

        elif not aln.aligned:
            unaligned_len.append(len(aln.read.seq))

    strandness = float(pos_strand) / num_aligned
    out_primary.close()
    unaligned_len = numpy.array(unaligned_len)
    return unaligned_len, strandness


def primary_and_unaligned_circular(sam_alnm_file, prefix, meta_list=None, ref_edge_max_dist=400, query_min_aln_len=100,
                                   include_other_primary_alns=True):
    '''
    Function to extract alignments of reads at extremities of circular genome references
    :param sam_alnm_file: path of input SAM file
    :param prefix: prefix of output SAM file
    :param meta_list: a dictionary of metagenome information, including circularity
    :param ref_edge_max_dist: max distance of alignments from extremities of references
    :param query_min_aln_len: min length of split alignments
    :param include_other_primary_alns: include other primary alignments
    :return: Outputs a sam file that contains: 1) primary alignments for all reads that can be aligned and 2)
             split alignments for reads that are aligned to both ends of a circular genome
             Return unaligned_len list, and strandness information
    '''

    # TODO: check meta_list, handle the case when a metagenome contains both circular and linear genomes

    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=False)
    tmp = open("tmp", 'w')
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0

    # extract the lengths of all reference sequences
    ref_lengths = dict()
    for info in in_sam_file.header['SQ']:
        ref_lengths[info['SN']] = info['LN']

    # info for the current read
    current_qname = None
    current_edge_alns = {}
    primary_edge = [False, False]
    current_primary_aln = None

    # parse alignment file, assuming that all alignments of a read are grouped together
    for aln in in_sam_file.fetch(until_eof=True):
        if aln.is_unmapped:
            unaligned_len.append(aln.query_length)
            continue
        if aln.query_name != current_qname:
            # a new read; evaluate all alignments for the previous read
            is_circular = False
            if current_qname == "ERR3152366.2599":
                print(current_edge_alns, primary_edge)
            if len(current_edge_alns) > 1 and any([x[0] for x in current_edge_alns.values()]) and \
                    any(x[1] for x in current_edge_alns.values()) and any(primary_edge):
                for a in current_edge_alns:
                    if a != current_primary_aln and (current_edge_alns[a][0] ^ primary_edge[0]) and \
                            (a.flag == aln_flag + 2048 or a.flag == aln_flag + 256):
                        is_circular = True
                        tmp.write(current_primary_aln.qname + '\t' + current_primary_aln.reference_name + '\t' +
                                  str(current_primary_aln.reference_start) + '\t' +
                                  str(current_primary_aln.reference_end) + '\t' + str(current_primary_aln.query_alignment_length) + '\t' +
                                  str(current_primary_aln.flag) + '\t' + str(current_primary_aln.is_secondary) + '\t' +
                                  str(current_primary_aln.is_supplementary) + '\t' +
                                  ','.join(str(x) for x in current_primary_aln.cigartuples[0]) + '\t' +
                                  ','.join(str(x) for x in current_primary_aln.cigartuples[-1]) + '\n')
                        out_sam_file.write(current_primary_aln)
                        tmp.write(a.qname + '\t' + a.reference_name + '\t' + str(a.reference_start) + '\t' +
                                  str(a.reference_end) + '\t' + str(a.query_alignment_length) + '\t' +
                                  str(a.flag) + '\t' + str(a.is_secondary) + '\t' +
                                  str(a.is_supplementary) + '\t' +
                                  ','.join(str(x) for x in a.cigartuples[0]) + '\t' +
                                  ','.join(str(x) for x in a.cigartuples[-1]) + '\n')
                        out_sam_file.write(a)
                        num_aligned += 1
                        if aln_flag == 0:
                            pos_strand += 1

                    else:
                        continue

            if include_other_primary_alns and not is_circular and current_primary_aln is not None:
                out_sam_file.write(current_primary_aln)
                num_aligned += 1
                if current_primary_aln.flag == 0:
                    pos_strand += 1

            # reset info for the new read
            current_qname = aln.query_name
            current_edge_alns = {}
            current_primary_aln = None
            primary_edge = [False, False]

        if not aln.is_secondary and not aln.is_supplementary:
            # a primary alignment
            current_primary_aln = aln
            aln_flag = aln.flag

        if aln.query_alignment_length >= query_min_aln_len and aln.reference_name == current_primary_aln.reference_name:
            if aln.reference_end >= ref_lengths[aln.reference_name] - 1 - ref_edge_max_dist:
                # read was aligned to reference end
                current_edge_alns[aln] = [False, True]
                if aln == current_primary_aln:
                    primary_edge[1] = True
            elif aln.reference_start <= ref_edge_max_dist:
                # read was aligned to reference start
                current_edge_alns[aln] = [True, False]
                if aln == current_primary_aln:
                    primary_edge[0] = True

    in_sam_file.close()

    # evaluate the final read
    is_circular = False
    if len(current_edge_alns) > 1 and any([x[0] for x in current_edge_alns.values()]) and \
            any(x[1] for x in current_edge_alns.values()) and any(primary_edge):
        for a in current_edge_alns:
            if a != current_primary_aln and (current_edge_alns[a][0] ^ primary_edge[0]) and \
                    (a.flag == aln_flag + 2048 or a.flag == aln_flag + 256):
                is_circular = True
                tmp.write(current_primary_aln.qname + '\t' + current_primary_aln.reference_name + '\t' +
                          str(current_primary_aln.reference_start) + '\t' +
                          str(current_primary_aln.reference_end) + '\t' + str(
                    current_primary_aln.query_alignment_length) + '\t' +
                          str(current_primary_aln.flag) + '\t' + str(current_primary_aln.is_secondary) + '\t' +
                          str(current_primary_aln.is_supplementary) + '\t' +
                          ','.join(str(x) for x in current_primary_aln.cigartuples[0]) + '\t' +
                          ','.join(str(x) for x in current_primary_aln.cigartuples[-1]) + '\n')
                out_sam_file.write(current_primary_aln)
                tmp.write(a.qname + '\t' + a.reference_name + '\t' + str(a.reference_start) + '\t' +
                          str(a.reference_end) + '\t' + str(a.query_alignment_length) + '\t' +
                          str(a.flag) + '\t' + str(a.is_secondary) + '\t' +
                          str(a.is_supplementary) + '\t' +
                          ','.join(str(x) for x in a.cigartuples[0]) + '\t' +
                          ','.join(str(x) for x in a.cigartuples[-1]) + '\n')
                out_sam_file.write(a)
                num_aligned += 1
                if aln_flag == 0:
                    pos_strand += 1

            else:
                continue

    if include_other_primary_alns and not is_circular and current_primary_aln is not None:
        out_sam_file.write(current_primary_aln)
        num_aligned += 1
        if aln.flag == 0:
            pos_strand += 1

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    out_sam_file.close()
    tmp.close()

    return unaligned_len, strandness
