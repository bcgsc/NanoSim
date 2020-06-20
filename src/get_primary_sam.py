#!/usr/bin/env python

from __future__ import with_statement
import HTSeq
import numpy
import pysam


def primary_and_unaligned(sam_alnm_file, prefix):
    out_primary = open(prefix + "_primary.sam", 'w')
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
    samfile = pysam.AlignmentFile(sam_alnm_file, "r")  # Read sam file using pysam for extracting sam header
    out_primary.write(str(samfile.header))  # Write sam header into out_primary
    sam_reader = HTSeq.SAM_Reader
    alignments = sam_reader(sam_alnm_file)
    for alnm in alignments:
        if alnm.aligned and (not alnm.not_primary_alignment and not alnm.supplementary):
            out_primary.write(alnm.get_sam_line() + '\n')  # Use get_sam_line() instead of original_sam_line, and add a newline)
            num_aligned += 1
            if alnm.flag == 0:
                pos_strand += 1

        elif not alnm.aligned:
            unaligned_len.append(len(alnm.read.seq))

    strandness = float(pos_strand) / num_aligned
    out_primary.close()
    unaligned_len = numpy.array(unaligned_len)
    return unaligned_len, strandness

