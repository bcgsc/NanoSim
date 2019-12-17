#!/usr/bin/env python

from __future__ import with_statement
import HTSeq
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
