#!/usr/bin/env python


import sys
import getopt
from time import strftime
import re


def main():
    infile = sys.argv[1]
    # print strftime("%Y-%m-%d %H:%M:%S") + ": Start!"

    with open(infile, 'r') as f:
        for line in f:
            new = line.strip().split()
            ref = new[6]
            new_line = next(f)
            new = new_line.strip().split()
            query = new[6]

            i = 0
            while i <= len(ref):
                if i == 0:
                    r_prev = ref[i]
                    i += 1
                    last_pos = [0, 1]
                    continue
                if i == len(ref) or (ref[i] != '-' and ref[i] != r_prev):
                    j = last_pos[1] - 1
                    while j >= 0 and query[j] == r_prev:
                        j -= 1
                    if j < last_pos[0]:
                        last_pos[0] = j + 1

                    j = last_pos[1]
                    while j < len(ref) and query[j] == r_prev:
                        j += 1
                    if j > last_pos[1]:
                        last_pos[1] = j

                    if last_pos[1] - last_pos[0] >= 6:
                        r_seq = ''.join(ref[last_pos[0]:last_pos[1]].split('-'))
                        r_homo = re.findall(r_prev + '+', r_seq)
                        r_len = max(len(h) for h in r_homo)

                        q_seq = ''.join(query[last_pos[0]:last_pos[1]].split('-'))
                        q_homo = re.findall(r_prev + '+', q_seq)
                        q_len = max(len(h) for h in q_homo) if len(q_homo) != 0 else 0

                        print ref[last_pos[0]:last_pos[1]], query[last_pos[0]:last_pos[1]], r_len, q_len

                    if j == len(ref):
                        break

                    i = last_pos[1]
                    while ref[i] == '-':
                        i += 1
                    last_pos[0] = i
                    r_prev = ref[i]
                    i += 1
                    last_pos[1] = i

                elif ref[i] == r_prev:
                    last_pos[1] += 1
                    i += 1
                elif ref[i] == '-':
                    i += 1

    # print strftime("%Y-%m-%d %H:%M:%S") + ": Finished!"


if __name__ == "__main__":
    main()

