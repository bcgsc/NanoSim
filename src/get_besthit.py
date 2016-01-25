#!/usr/bin/env python

from __future__ import with_statement


def besthit_and_unaligned(infile, outmaf, outfile):
    align_dict = {}
    out1 = open(outfile + "_besthit.maf", 'w')
    unaligned_dict = {}

    with open(outmaf, 'r') as f:
        for line in f:
            new_line = next(f)
            new_info = new_line.strip().split()
            if new_info[1] not in align_dict:
                align_dict[new_info[1]] = [int(new_info[3]), new_line, False]
            else:
                if align_dict[new_info[1]][0] < int(new_info[3]):
                    align_dict[new_info[1]] = [int(new_info[3]), new_line, False]

    with open(outmaf, 'r') as f1:
        for line in f1:
            ref = line
            new_line = next(f1)
            new_info = new_line.split()
            name = new_info[1]
            length = int(new_info[3])
            if align_dict[name][0] == length and not align_dict[name][2]:
                out1.write(ref + new_line)
                align_dict[name][2] = True

    with open(infile, 'r') as f2:
        for line in f2:
            if line[0] == ">":
                name = line.strip()[1:]
                flag = False
                if name not in align_dict:
                    last_name = name
                    flag = True
            else:
                if flag:
                    unaligned_dict[last_name] = len(line.strip())

    out1.close()
    return unaligned_dict.values()
