#!/usr/bin/env python

from __future__ import with_statement
import HTSeq


def primary_and_unaligned(sam_alnm_file, outfile):
    out_sam = outfile + ".sam"
    if out_sam == sam_alnm_file:
        out_sam = outfile + "_processed.sam"
    out1 = open(out_sam, 'w')
    out2 = open(outfile + "_primary.sam", 'w')
    unaligned_dict = []
    sam_reader = HTSeq.SAM_Reader
    alignments = sam_reader(sam_alnm_file)
    for alnm in alignments:
        out1.write(alnm.original_sam_line)
        if alnm.aligned and not alnm.not_primary_alignment and not alnm.supplementary:
            out2.write(alnm.original_sam_line)
        else:
            unaligned_dict.append(len(alnm.read.seq))

    out1.close()
    out2.close()
    return unaligned_dict
