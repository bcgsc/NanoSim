#!/usr/bin/env python

from __future__ import with_statement
import numpy


def add_dict(error, dic):
    if error not in dic:
        last_elememt = len(dic)
        for i in xrange(last_elememt, error + 1):
            dic[i] = 0
    dic[error] += 1


def add_match(prev, succ, match_list):
    # expand the match_list matrix to the biggest possible size
    expand = max(prev, succ) + 1
    if expand > len(match_list):
        last_element = len(match_list)
        for i in xrange(0, last_element):
            for j in xrange(last_element, expand):
                match_list[i][j] = 0
        for i in xrange(last_element, expand):
            match_list[i] = {}
            for j in xrange(0, expand):
                match_list[i][j] = 0

    match_list[prev][succ] += 1


def hist(outfile):
    out_match = open(outfile + "_match.hist", 'w')
    out_mis = open(outfile + "_mis.hist", 'w')
    out_ins = open(outfile + "_ins.hist", 'w')
    out_del = open(outfile + "_del.hist", 'w')
    out1 = open(outfile + "_error_markov_model", 'w')
    out2 = open(outfile + "_match_markov_model", 'w')
    out3 = open(outfile + "_first_match.hist", 'w')

    dic_match = {}
    dic_first_match = {}
    dic_mis = {}
    dic_ins = {}
    dic_del = {}
    match_list = {}
    match_bin = {}
    error_list = {"mis/mis": 0, "mis/ins": 0, "mis/del": 0, "ins/mis": 0, "ins/ins": 0, "ins/del": 0,
                  "del/mis": 0, "del/ins": 0, "del/del": 0, "mis0/mis": 0, "mis0/ins": 0, "mis0/del": 0,
                  "del0/mis": 0, "del0/ins": 0, "del0/del": 0, "ins0/mis": 0, "ins0/ins": 0, "ins0/del": 0}
    first_error = {"mis": 0, "ins": 0, "del": 0}

    for x in xrange(0, 150):
        dic_match[x] = 0
        match_list[x] = {}
        for y in xrange(0, 150):
            match_list[x][y] = 0
        for key in match_bin.keys():
            match_bin[key][x] = 0
    for x in xrange(0, 150):
        dic_first_match[x] = 0

    for x in xrange(0, 30):
        dic_mis[x] = 0
        dic_ins[x] = 0
        dic_del[x] = 0

    with open(outfile + "_besthit.maf", 'r') as f:
        for line in f:
            prev_match = 0
            prev_error = ""
            flag = True
            new = line.strip().split()
            ref = new[6]
            new_line = next(f)
            new = new_line.strip().split()
            query = new[6]
            match = 0
            mismatch = 0
            ins = 0
            dele = 0
            for i in xrange(0, len(ref)):
                if ref[i] == query[i]:
                    if mismatch != 0:
                        add_dict(mismatch, dic_mis)
                        mismatch = 0
                        if flag:
                            flag = False
                            first_error["mis"] += 1
                        else:
                            error_list[prev_error + "/" + "mis"] += 1
                        prev_error = "mis"
                    elif ins != 0:
                        add_dict(ins, dic_ins)
                        ins = 0
                        if flag:
                            flag = False
                            first_error["ins"] += 1
                        else:
                            error_list[prev_error + "/" + "ins"] += 1
                        prev_error = "ins"
                    elif dele != 0:
                        add_dict(dele, dic_del)
                        dele = 0
                        if flag:
                            flag = False
                            first_error["del"] += 1
                        else:
                            error_list[prev_error + "/" + "del"] += 1
                        prev_error = "del"
                    match += 1
                    if i == len(ref) - 1 and match != 0:
                        add_match(prev_match, match, match_list)
                elif ref[i] == '-':
                    if match != 0:
                        if flag:
                            add_dict(match, dic_first_match)
                            prev_match = match
                        else:
                            add_dict(match, dic_match)
                            add_match(prev_match, match, match_list)
                            prev_match = match
                        match = 0
                    elif mismatch != 0:
                        add_dict(mismatch, dic_mis)
                        dic_match[0] += 1
                        add_match(prev_match, 0, match_list)
                        prev_match = 0
                        mismatch = 0
                        if flag:
                            flag = False
                            first_error["mis"] += 1
                        else:
                            error_list[prev_error + "/" + "mis"] += 1
                        prev_error = "mis0"
                    ins += 1
                elif query[i] == '-':
                    if match != 0:
                        if flag:
                            add_dict(match, dic_first_match)
                            prev_match = match
                        else:
                            add_dict(match, dic_match)
                            add_match(prev_match, match, match_list)
                            prev_match = match
                        match = 0
                    elif mismatch != 0:
                        add_dict(mismatch, dic_mis)
                        dic_match[0] += 1
                        add_match(prev_match, 0, match_list)
                        prev_match = 0
                        mismatch = 0
                        if flag:
                            flag = False
                            first_error["mis"] += 1
                        else:
                            error_list[prev_error + "/" + "mis"] += 1
                        prev_error = "mis0"
                    dele += 1
                else:
                    if match != 0:
                        if flag:
                            add_dict(match, dic_first_match)
                            prev_match = match
                        else:
                            add_dict(match, dic_match)
                            add_match(prev_match, match, match_list)
                            prev_match = match
                        match = 0
                    elif ins != 0:
                        add_dict(ins, dic_ins)
                        add_dict(match, dic_match)
                        add_match(prev_match, 0, match_list)
                        prev_match = 0
                        ins = 0
                        if flag:
                            flag = False
                            first_error["ins"] += 1
                        else:
                            error_list[prev_error + "/" + "ins"] += 1
                        prev_error = "ins0"
                    elif dele != 0:
                        add_dict(dele, dic_del)
                        add_dict(match, dic_match)
                        add_match(prev_match, 0, match_list)
                        prev_match = 0
                        dele = 0
                        if flag:
                            flag = False
                            first_error["del"] += 1
                        else:
                            error_list[prev_error + "/" + "del"] += 1
                        prev_error = "del0"
                    mismatch += 1

    # write the histogram for other matches and errors:
    out_match.write("number of bases\tMatches:\n")
    for key in dic_match:
        out_match.write(str(key) + "\t" + str(dic_match[key]) + "\n")
    out_match.close()

    out_mis.write("number of bases\tMismatches:\n")
    for key in dic_mis:
        out_mis.write(str(key) + "\t" + str(dic_mis[key]) + "\n")
    out_mis.close()

    out_ins.write("number of bases\tInsertions:\n")
    for key in dic_ins:
        out_ins.write(str(key) + "\t" + str(dic_ins[key]) + "\n")
    out_ins.close()

    out_del.write("number of bases\tDeletions:\n")
    for key in dic_del:
        out_del.write(str(key) + "\t" + str(dic_del[key]) + "\n")
    out_del.close()

    predecessor = {"mis": error_list["mis/mis"] + error_list["mis/ins"] + error_list["mis/del"],
                   "ins": error_list["ins/mis"] + error_list["ins/ins"] + error_list["ins/del"],
                   "del": error_list["del/mis"] + error_list["del/ins"] + error_list["del/del"],
                   "mis0": error_list["mis0/mis"] + error_list["mis0/ins"] + error_list["mis0/del"],
                   "ins0": error_list["ins0/mis"] + error_list["ins0/ins"] + error_list["ins0/del"],
                   "del0": error_list["del0/mis"] + error_list["del0/ins"] + error_list["del0/del"]}
    out1.write("succedent \tmis\tins\tdel\n")

    num_of_first = sum(first_error.values())
    out1.write("start\t" + str(first_error["mis"] * 1.0 / num_of_first) + "\t" + str(first_error["ins"] * 1.0 / num_of_first) +
               "\t" + str(first_error["del"] * 1.0 / num_of_first))
    for x in ["mis", "ins", "del", "mis0", "ins0", "del0"]:
        out1.write("\n" + x)
        for y in ["mis", "ins", "del"]:
            if predecessor[x] == 0:
                out1.write("\t" + "0")
            else:
                out1.write("\t" + str(error_list[x + "/" + y] * 1.0 / predecessor[x]))

    # Match markov model
    count = 0
    for k1 in sorted(match_list.keys()):
        count += sum(match_list[k1].values())
    # 15 bins for the precedent in the match pair, each bin has roughly count/15 match events
    bin_size = count / 15
    k_of_bin = 0
    k_of_match_list = 0
    last_k = 0
    count_each_bin = {}
    while k_of_bin < 15:
        if k_of_match_list >= len(match_list):
            break
        match_bin[k_of_bin] = {}
        for i in xrange(0, len(match_list)):
            match_bin[k_of_bin][i] = 0
        tmp_count = 0
        while tmp_count < bin_size and k_of_match_list < len(match_list):
            new_added = sum(match_list[k_of_match_list].values())
            if abs(tmp_count + new_added - bin_size) > abs(tmp_count - bin_size) and tmp_count != 0:
                break
            else:
                tmp_count += new_added
                k_of_match_list += 1
        for k1 in xrange(last_k, k_of_match_list):
            for k2 in xrange(0, len(match_list[k1])):
                match_bin[k_of_bin][k2] += match_list[k1][k2]

        count_each_bin[k_of_bin] = [(last_k, k_of_match_list), tmp_count]
        last_k = k_of_match_list
        k_of_bin += 1

    if k_of_match_list < len(match_list):
        tmp_count = 0
        for k1 in xrange(last_k, len(match_list)):
            tmp_count += sum(match_list[k1].values())
            for k2 in xrange(0, len(match_list)):
                match_bin[k_of_bin - 1][k2] += match_list[k1][k2]
        count_each_bin[k_of_bin - 1][1] += tmp_count

    count_prob = [0] * len(match_bin)
    out2.write("bins\t" + "\t".join("%s-%s" % tup[0] for tup in count_each_bin.values()) + '\n')
    for i in xrange(0, len(match_list)):
        out2.write(str(i) + "-" + str(i + 1))
        for k_of_bin in match_bin:
            if count_each_bin[k_of_bin][1] == 0:
                out2.write("\t" + "0")
            else:
                count_prob[k_of_bin] += match_bin[k_of_bin][i] * 1.0 / count_each_bin[k_of_bin][1]
                out2.write("\t" + str(count_prob[k_of_bin]))
        out2.write('\n')

    out2.close()

    # First match profile:
    out3.write("bin\t0-50000\n")
    count_prob = 0
    total_first_match = sum(dic_first_match.values())
    for i in xrange(0, len(dic_first_match)):
        count_prob += dic_first_match[i] * 1.0 / total_first_match
        out3.write(str(i) + "-" + str(i + 1) + "\t" + str(count_prob) + '\n')

    out3.close()
