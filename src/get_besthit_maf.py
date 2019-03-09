#!/usr/bin/env python

from __future__ import with_statement
import numpy


def besthit_and_unaligned(infile, outmaf, prefix):
    align_dict = {}
    out_best = open(prefix + "_besthit.maf", 'w')
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
    with open(outmaf, 'r') as f:
        for line in f:
            query = next(f)
            query_info = query.strip().split()
            if query_info[1] not in align_dict:
                align_dict[query_info[1]] = [int(query_info[3]), query, False]
            else:
                if align_dict[query_info[1]][0] < int(query_info[3]):
                    align_dict[query_info[1]] = [int(query_info[3]), query, False]

    with open(outmaf, 'r') as f1:
        for line in f1:
            ref = line
            ref_info = ref.split()
            ref_strand = ref_info[4]
            query = next(f1)
            query_info = query.split()
            name = query_info[1]
            query_strand = query_info[4]
            length = int(query_info[3])
            if align_dict[name][0] == length and not align_dict[name][2]:
                out_best.write(ref + query)
                align_dict[name][2] = True
                num_aligned += 1
                if ref_strand == query_strand:
                    pos_strand += 1

    strandness = float(pos_strand) / num_aligned

    with open(infile, 'r') as f2:
        for line in f2:
            if line[0] == ">":
                name = line.strip().split()[0][1:]
                flag = False
                if name not in align_dict:
                    last_name = name
                    flag = True
            else:
                if flag:
                    unaligned_len.append(len(line.strip()))

    out_best.close()
    unaligned_len = numpy.array(unaligned_len)
    return unaligned_len, strandness
