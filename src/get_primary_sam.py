#!/usr/bin/env python

from __future__ import with_statement
import HTSeq


def primary_and_unaligned(outsam, outfile):
    out1 = open(outfile + "_primary.sam", 'w')
    unaligned_dict = []
    sam_reader = HTSeq.SAM_Reader
    alignments = sam_reader(outsam)
    for alnm in alignments:
        if alnm.aligned and not alnm.not_primary_alignment and not alnm.supplementary:
            out1.write(alnm.original_sam_line)
        else:
            unaligned_dict.append(len(alnm.read.seq))

    out1.close()
    return unaligned_dict
