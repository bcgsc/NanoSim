#!/usr/bin/env python

from __future__ import with_statement
import HTSeq


def primary_and_unaligned(sam_alnm_file, prefix):
    out_primary = open(prefix + "_primary.sam", 'w')
    unaligned_len = []
    sam_reader = HTSeq.SAM_Reader
    alignments = sam_reader(sam_alnm_file)
    for alnm in alignments:
        if alnm.aligned and (not alnm.not_primary_alignment and not alnm.supplementary):
            out_primary.write(alnm.original_sam_line)
        elif not alnm.aligned:
            unaligned_len.append(len(alnm.read.seq))

    out_primary.close()
    return unaligned_len
