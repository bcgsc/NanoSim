#!/usr/bin/env python

"""
Written by: Ka Ming Nip (BC Cancer Genome Sciences Centre)
Adapted by: Chen Yang (BC Cancer Genome Sciences Centre)
"""

from __future__ import with_statement
import sys
import pysam
import numpy


def primary_and_unaligned(sam_alnm_file, prefix, meta_list, ref_edge_max_dist=400, query_min_aln_len=100,
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
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
        
    # extract the lengths of all reference sequences
    ref_lengths = dict()
    for info in in_sam_file.header['SQ']:
        ref_lengths[info['SN']] = info['LN']

    # info for the current read
    current_qname = None
    current_start_rnames = set()
    current_end_rnames = set()
    current_edge_alns = []
    current_primary_aln = None

    # parse alignment file, assuming that all alignments of a read are grouped together
    for aln in in_sam_file.fetch(until_eof=True):
        if not aln.is_unmapped:
            if aln.query_name != current_qname:
                # a new read; evaluate all alignments for the previous read
                is_circular = False
                if len(current_start_rnames) > 0 and len(current_end_rnames) > 0:
                    for rname in (current_start_rnames & current_end_rnames):
                        is_circular = True
                        # write alignments for this reference
                        for a in current_edge_alns:
                            flag = False
                            if a.reference_name == rname:
                                out_sam_file.write(a)
                                if not flag:
                                    num_aligned += 1
                                    if aln.flag == 0:
                                        pos_strand += 1
                                    flag = True
                
                if include_other_primary_alns and not is_circular and current_primary_aln is not None:
                    out_sam_file.write(current_primary_aln)
                    num_aligned += 1
                    if aln.flag == 0:
                        pos_strand += 1
                        
                # reset info for the new read
                current_qname = aln.query_name
                current_start_rnames = set()
                current_end_rnames = set()
                current_edge_alns = []
                current_primary_aln = None

            if not aln.is_secondary and not aln.is_supplementary:
                # a primary alignment
                current_primary_aln = aln
            
            if aln.query_alignment_length >= query_min_aln_len:
                if aln.reference_end >= ref_lengths[aln.reference_name] - 1 - ref_edge_max_dist:
                    # read was aligned to reference end
                    current_end_rnames.add(aln.reference_name)
                    current_edge_alns.append(aln)
                elif aln.reference_start <= ref_edge_max_dist:
                    # read was aligned to reference start
                    current_start_rnames.add(aln.reference_name)
                    current_edge_alns.append(aln)

        else:
            unaligned_len.append(aln.query_length)
    
    in_sam_file.close()
    
    # evaluate the final read
    is_circular = False
    if len(current_start_rnames) > 0 and len(current_end_rnames) > 0:
        for rname in (current_start_rnames & current_end_rnames):
            is_circular = True
            # write alignments for this reference
            for a in current_edge_alns:
                flag = False
                if a.reference_name == rname:
                    out_sam_file.write(aln)
                    if not flag:
                        num_aligned += 1
                        if aln.flag == 0:
                            pos_strand += 1
                        flag = True
    
    if include_other_primary_alns and not is_circular and current_primary_aln is not None:
        out_sam_file.write(current_primary_aln)
        num_aligned += 1
        if aln.flag == 0:
            pos_strand += 1

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned
    
    out_sam_file.close()

    return unaligned_len, strandness


