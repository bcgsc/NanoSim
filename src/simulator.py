#!/usr/bin/env python
"""
@author: Chen Yang & Saber HafezQorani
This script generates simulated Oxford Nanopore 2D reads (genomic and transcriptomic - cDNA/directRNA).
"""

from __future__ import print_function
from __future__ import with_statement

import multiprocessing as mp
from subprocess import call
from textwrap import dedent
import sys
import os
import HTSeq
import pysam
import numpy
import random
import re
import copy
import argparse
import joblib
from time import strftime
import numpy as np

if sys.version_info[0] < 3:
    from string import maketrans
    trantab = maketrans("T", "U")
else:
    trantab = str.maketrans("T", "U")

try:
    from six.moves import xrange
except ImportError:
    pass
import mixed_model as mm
import norm_distr as nd

PYTHON_VERSION = sys.version_info
VERSION = "2.5.0"
PRORAM = "NanoSim"
AUTHOR = "Chen Yang, Saber Hafezqorani (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca; shafezqorani@bcgsc.ca"

BASES = ['A', 'T', 'C', 'G']


def select_ref_transcript(input_dict):
    length = 0
    while True:
        p = random.random()
        for key, val in input_dict.items():
            if key[0] <= p < key[1]:
                length = val[1]
                break
        if length != 0:
            break
    return val[0], length


def list_to_range(input_list, min_l):
    l = [min_l]
    l.extend(input_list)
    output_list = []
    for i in xrange(0, len(l) - 1):
        r = (l[i], l[i+1])
        output_list.append(r)
    return output_list


def make_cdf(dict_exp, dict_len):
    sum_exp = 0
    list_value = []
    for item in dict_exp:
        if item in dict_len:
            sum_exp += dict_exp[item]
    for item in dict_exp:
        if item in dict_len:
            value = dict_exp[item] / float(sum_exp)
            list_value.append((item, value))

    sorted_value_list = sorted(list_value, key=lambda x: x[1])
    sorted_only_values = [x[1] for x in sorted_value_list]
    list_cdf = numpy.cumsum(sorted_only_values)
    ranged_cdf_list = list_to_range(list_cdf, 0)

    ecdf_dict = {}
    for i in xrange(len(ranged_cdf_list)):
        cdf_range = ranged_cdf_list[i]
        ecdf_dict[cdf_range] = (sorted_value_list[i][0], dict_len[sorted_value_list[i][0]])
        # ecdf_dict[cdf_range] = dict_len[sorted_value_list[i][0]]

    return ecdf_dict


def ref_len_from_structure(input):
    l = 0
    for item in input:
        if item[0] == "exon":
            l += item[-1]
    return l


def select_nearest_kde2d(sampled_2d_lengths, ref_len_total):
    fc = sampled_2d_lengths[:, 0]
    idx = np.abs(fc - ref_len_total).argmin()
    return int(sampled_2d_lengths[idx][1])


def update_structure(ref_trx_structure, IR_markov_model):
    count = 0
    for item in ref_trx_structure:
        if item[0] == "intron":
            count += 1

    list_states = []
    flag_ir = False
    prev_state = "start"
    for i in range (0, count):
        p = random.random()
        for key in IR_markov_model[prev_state]:
            if key[0] <= p < key[1]:
                flag = IR_markov_model[prev_state][key]
                if flag == "IR":
                    flag_ir = True
                list_states.append(flag)
                prev_state = flag

    if flag_ir:
        ref_trx_structure_temp = copy.deepcopy(ref_trx_structure)
        j = -1
        for i in xrange (0, len(ref_trx_structure_temp)):
            if ref_trx_structure_temp[i][0] == "intron":
                j += 1
                if list_states[j] == "IR":
                    ref_trx_structure_temp[i] = ("retained_intron",) + ref_trx_structure_temp[i][1:]
    else:
        ref_trx_structure_temp = ref_trx_structure

    return flag_ir, ref_trx_structure_temp


def extract_read_pos(length, ref_len, ref_trx_structure):
    # The aim is to create a genomic interval object
    # example: iv = HTSeq.GenomicInterval( "chr3", 123203, 127245, "+" )
    chrom = ""
    start = 0
    end = 0
    strand = ""
    list_intervals = []
    start_pos = random.randint(0, ref_len - length)
    flag = False
    for item in ref_trx_structure:
        chrom = item[1]
        if item[0] == "exon":
            if not flag:  # if it is first exon that I started extracting from
                if start_pos < item[-1]:
                    # print ("1", item, start_pos, length, start, end)
                    flag = True
                    start = start_pos + item[2]
                    if (start + length) < item[3]:
                        end = start + length
                        length = 0
                    else:
                        end = item[3]
                        length -= (item[3] - start)
                    start_pos = 0
                else:
                    # print("2", item, start_pos, length, start, end)
                    start_pos -= item[-1]
            else:  # if it is NOT the first exon that I start extracting from
                # print("3", item, start_pos, length, start, end)
                start_pos = 0
                if (item[2] + length) < item[3]:
                    end = item[2] + length
                    length = 0

                else:
                    end = item[3]
                    length -= (item[-1])

        elif item[0] == "retained_intron":
            # print("4", item, start_pos, length, start, end)
            if flag != False:
                end = item[3]
        elif item[0] == "intron":
            # print("5", item, start_pos, length, start, end)
            iv = HTSeq.GenomicInterval(chrom, start, end, ".")
            if iv.length != 0:
                list_intervals.append(iv)
            chrom = ""
            start = 0
            end = 0
            strand = ""
            start_pos = 0
            flag = False
    iv = HTSeq.GenomicInterval(chrom, start, end, ".")
    if iv.length != 0:
        list_intervals.append(iv)

    return list_intervals


def read_ecdf(profile):
    # We need to count the number of zeros. If it's over 10 zeros, l_len/l_ratio need to be changed to higher.
    # Because it's almost impossible that the ratio is much lower than the lowest historical value.
    header = profile.readline()
    header_info = header.strip().split()
    ecdf_dict = {}
    lanes = len(header_info[1:])

    for i in header_info[1:]:
        boundaries = i.split('-')
        ecdf_dict[(int(boundaries[0])), int(boundaries[1])] = {}

    ecdf_key = sorted(ecdf_dict.keys())
    l_prob = [0.0] * lanes
    l_ratio = [0.0] * lanes

    for line in profile:
        new = line.strip().split('\t')
        ratio = [float(x) for x in new[0].split('-')]
        prob = [float(x) for x in new[1:]]
        for i in xrange(lanes):
            if prob[i] == l_prob[i]:
                continue
            else:
                if l_prob[i] != 0:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] = (l_ratio[i], ratio[1])
                else:
                    ecdf_dict[ecdf_key[i]][(l_prob[i], prob[i])] \
                        = (max(l_ratio[i], ratio[1] - 10 * (ratio[1] - ratio[0])), ratio[1])
                l_ratio[i] = ratio[1]
                l_prob[i] = prob[i]

    for i in xrange(0, len(ecdf_key)):
        last_key = sorted(ecdf_dict[ecdf_key[i]].keys())[-1]
        last_value = ecdf_dict[ecdf_key[i]][last_key]
        ecdf_dict[ecdf_key[i]][last_key] = (last_value[0], ratio[1])

    return ecdf_dict


def get_length_kde(kde, num, log=False, flatten=True):
    tmp_list = kde.sample(num)
    if log:
        tmp_list = np.power(10, tmp_list) - 1
    length_list = tmp_list.flatten()
    if not flatten:
        return tmp_list
    else:
        return length_list


def read_profile(ref_g, ref_t, number, model_prefix, per, mode, strandness, exp, model_ir, dna_type):
    global number_aligned, number_unaligned
    global match_ht_list, error_par, trans_error_pr, match_markov_model
    global kde_aligned, kde_ht, kde_ht_ratio, kde_unaligned, kde_aligned_2d
    global seq_dict, seq_len, max_chrom
    global strandness_rate

    if mode == "genome":
        global genome_len
        ref = ref_g
    else:
        global dict_ref_structure, dict_exp, ecdf_dict_ref_exp, genome_fai
        ref = ref_t

    if strandness is None:
        with open(model_prefix + "_strandness_rate", 'r') as strand_profile:
            strandness_rate = float(strand_profile.readline().split("\t")[1])
    else:
        strandness_rate = strandness

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference \n")
    sys.stdout.flush()
    seq_dict = {}
    seq_len = {}

    # Read in the reference genome/transcriptome
    with open(ref, 'r') as infile:
        for seqN, seqS, seqQ in readfq(infile):
            info = re.split(r'[_\s]\s*', seqN)
            chr_name = "-".join(info)
            seq_dict[chr_name.split(".")[0]] = seqS
            seq_len[chr_name.split(".")[0]] = len(seqS)

    max_chrom = max(seq_len.values())
    if mode == "genome":
        genome_len = sum(seq_len.values())
        if len(seq_dict) > 1 and dna_type == "circular":
            sys.stderr.write("Do not choose circular if there is more than one chromosome in the genome!")
            sys.exit(1)
    else:
        if model_ir:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference genome and create .fai index file\n")
            sys.stdout.flush()
            # create and read the .fai file of the reference genome
            genome_fai = pysam.Fastafile(ref_g)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in expression profile\n")
        sys.stdout.flush()
        dict_exp = {}
        with open(exp, 'r') as exp_file:
            header = exp_file.readline()
            for line in exp_file:
                parts = line.split("\t")
                transcript_id = parts[0].split(".")[0]
                tpm = float(parts[2])
                if transcript_id.startswith("ENS") and tpm > 0:
                    dict_exp[transcript_id] = tpm

        # create the ecdf dict considering the expression profiles
        ecdf_dict_ref_exp = make_cdf(dict_exp, seq_len)

        if model_ir:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in IR markov model\n")
            sys.stdout.flush()

            global IR_markov_model
            IR_markov_model = {}
            with open(model_prefix + "_IR_markov_model", "r") as  IR_markov:
                IR_markov.readline()
                for line in IR_markov:
                    info = line.strip().split()
                    k = info[0]
                    IR_markov_model[k] = {}
                    IR_markov_model[k][(0, float(info[1]))] = "no_IR"
                    IR_markov_model[k][(float(info[1]), float(info[1]) + float(info[2]))] = "IR"

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in GFF3 annotation file\n")
            sys.stdout.flush()
            dict_ref_structure = {}
            gff_file = model_prefix + "_added_intron_final.gff3"
            gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)
            for feature in gff_features:
                if feature.type == "exon" or feature.type == "intron":
                    if "transcript_id" in feature.attr:
                        feature_id = feature.attr['transcript_id']
                    elif "Parent" in feature.attr:
                        info = feature.name.split(":")
                        if len(info) == 1:
                            feature_id = info[0]
                        else:
                            if info[0] == "transcript":
                                feature_id = info[1]
                            else:
                                continue
                    else:
                        continue

                    feature_id = feature_id.split(".")[0]
                    if feature_id not in dict_ref_structure:
                        dict_ref_structure[feature_id] = []

                    # remove "chr" from chromosome names to be constant
                    if "chr" in feature.iv.chrom:
                        feature.iv.chrom = feature.iv.chrom.strip("chr")

                    dict_ref_structure[feature_id].append(
                        (feature.type, feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.length))

    if per:  # if parameter perfect is used, all reads should be aligned, number_aligned equals total number of reads
        number_aligned = number
    else:
        # Read model profile for match, mismatch, insertion and deletions
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read error profile\n")
        sys.stdout.flush()

        error_par = {}
        model_profile = model_prefix + "_model_profile"
        with open(model_profile, 'r') as mod_profile:
            mod_profile.readline()
            for line in mod_profile:
                new_line = line.strip().split("\t")
                if "mismatch" in line:
                    error_par["mis"] = [float(x) for x in new_line[1:]]
                elif "insertion" in line:
                    error_par["ins"] = [float(x) for x in new_line[1:]]
                else:
                    error_par["del"] = [float(x) for x in new_line[1:]]

        trans_error_pr = {}
        with open(model_prefix + "_error_markov_model", "r") as error_markov:
            error_markov.readline()
            for line in error_markov:
                info = line.strip().split()
                k = info[0]
                trans_error_pr[k] = {}
                trans_error_pr[k][(0, float(info[1]))] = "mis"
                trans_error_pr[k][(float(info[1]), float(info[1]) + float(info[2]))] = "ins"
                trans_error_pr[k][(1 - float(info[3]), 1)] = "del"

        with open(model_prefix + "_first_match.hist", 'r') as fm_profile:
            match_ht_list = read_ecdf(fm_profile)

        with open(model_prefix + "_match_markov_model", 'r') as mm_profile:
            match_markov_model = read_ecdf(mm_profile)

        # Read length of unaligned reads
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read KDF of unaligned reads\n")
        sys.stdout.flush()

        with open(model_prefix + "_reads_alignment_rate", 'r') as u_profile:
            new = u_profile.readline().strip()
            rate = new.split('\t')[1]
            if rate == "100%":
                number_aligned = number
            else:
                number_aligned = int(round(number * float(rate) / (float(rate) + 1)))
            number_unaligned = number - number_aligned

        if number_unaligned > 0:
            kde_unaligned = joblib.load(model_prefix + "_unaligned_length.pkl")

    # Read profile of aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read KDF of aligned reads\n")
    sys.stdout.flush()

    # Read ht profile
    kde_ht = joblib.load(model_prefix + "_ht_length.pkl")

    # Read head/unaligned region ratio
    kde_ht_ratio = joblib.load(model_prefix + "_ht_ratio.pkl")

    # Read length of aligned reads
    # If "perfect" is chosen, just use the total length ecdf profile, else use the length of aligned region on reference
    if per:
        kde_aligned = joblib.load(model_prefix + "_aligned_reads.pkl")
    else:
        if mode == "genome":
            kde_aligned = joblib.load(model_prefix + "_aligned_region.pkl")
        else:
            kde_aligned_2d = joblib.load(model_prefix + "_aligned_region_2d.pkl")


def mutate_homo(seq, k, basecaller, read_type):
    hp_arr = []  # [[base, start, end], ...]
    hp_length_hist = {}  # {length: {A/T: count, C/G: count} ...}
    hp_samples = {}  # {length: {A/T: [sample], C/G: [sample]} ...}

    # Finding homopolymers in sequence
    pattern = "A{" + re.escape(str(k)) + ",}|C{" + re.escape(str(k)) + ",}|G{" + re.escape(
        str(k)) + ",}|T{" + re.escape(str(k)) + ",}"

    for match in re.finditer(pattern, seq):
        hp_start = match.start()
        hp_end = match.end()
        length = hp_end - hp_start
        base = match.group()[0]
        hp_arr.append([base, hp_start, hp_end])
        if length not in hp_length_hist.keys():
            hp_length_hist[length] = {"A": 0, "T": 0, "C": 0, "G": 0}

        if base == "A":
            key = "A"
        elif base == "T":
            key = "T"
        elif base == "C":
            key = "C"
        else:
            key = "G"

        hp_length_hist[length][key] += 1

    # Obtaining samples from normal distributions
    for length in hp_length_hist.keys():
        hp_samples[length] = {}
        a_mu, a_sigma, t_mu, t_sigma, c_mu, c_sigma, g_mu, g_sigma = nd.get_nd_par(length, read_type, basecaller)

        if hp_length_hist[length]["A"] > 0:
            hp_samples[length]["A"] = np.random.normal(a_mu, a_sigma, hp_length_hist[length]["A"])
        if hp_length_hist[length]["T"] > 0:
            hp_samples[length]["T"] = np.random.normal(t_mu, t_sigma, hp_length_hist[length]["T"])
        if hp_length_hist[length]["C"] > 0:
            hp_samples[length]["C"] = np.random.normal(c_mu, c_sigma, hp_length_hist[length]["C"])
        if hp_length_hist[length]["G"] > 0:
            hp_samples[length]["G"] = np.random.normal(g_mu, g_sigma, hp_length_hist[length]["G"])

    # Mutating homopolymers in given sequence
    last_pos = 0
    mutated_seq = ""
    total_hp_size_change = 0
    mis_rate = nd.get_hpmis_rate(read_type, basecaller)
    for hp_info in hp_arr:
        base = hp_info[0]
        ref_hp_start = hp_info[1]
        ref_hp_end = hp_info[2]

        if base == "A":
            size = round(hp_samples[ref_hp_end - ref_hp_start]["A"][-1])
            hp_samples[ref_hp_end - ref_hp_start]["A"] = hp_samples[ref_hp_end - ref_hp_start]["A"][:-1]
        elif base == "T":
            size = round(hp_samples[ref_hp_end - ref_hp_start]["T"][-1])
            hp_samples[ref_hp_end - ref_hp_start]["T"] = hp_samples[ref_hp_end - ref_hp_start]["T"][:-1]
        elif base == "C":
            size = round(hp_samples[ref_hp_end - ref_hp_start]["C"][-1])
            hp_samples[ref_hp_end - ref_hp_start]["C"] = hp_samples[ref_hp_end - ref_hp_start]["C"][:-1]
        else:
            size = round(hp_samples[ref_hp_end - ref_hp_start]["G"][-1])
            hp_samples[ref_hp_end - ref_hp_start]["G"] = hp_samples[ref_hp_end - ref_hp_start]["G"][:-1]

        mutated_hp = base * int(size)
        mutated_hp_with_mis = ""
        for old_base in mutated_hp:
            p = random.random()
            if 0 < p <= mis_rate:
                tmp_bases = list(BASES)
                new_base = old_base
                while new_base != old_base:
                    new_base = random.choice(tmp_bases)
                mutated_hp_with_mis += new_base
            else:
                mutated_hp_with_mis += old_base
        mutated_seq = mutated_seq + seq[last_pos: ref_hp_start] + mutated_hp

        total_hp_size_change += int(size) - (ref_hp_end - ref_hp_start)
        last_pos = ref_hp_end

    return mutated_seq + seq[last_pos:]


# Taken from https://github.com/lh3/readfq
def readfq(fp):  # this is a generator function
    last = None  # this is a buffer keeping the last unprocessed line
    while True:  # mimic closure; is it a bad idea?
        if not last:  # the first record or a record following a fastq
            for l in fp:  # search for the start of the next record
                if l[0] in '>@':  # fasta/q header line
                    last = l[:-1]  # save this line
                    break
        if not last:
            break
        name, seqs, last = last[1:].partition(" ")[0], [], None
        for l in fp:  # read the sequence
            if l[0] in '@+>':
                last = l[:-1]
                break
            seqs.append(l[:-1])
        if not last or last[0] != '+':  # this is a fasta record
            yield name, ''.join(seqs), None  # yield a fasta record
            if not last:
                break
        else:  # this is a fastq record
            seq, leng, seqs = ''.join(seqs), 0, []
            for l in fp:  # read the quality
                seqs.append(l[:-1])
                leng += len(l) - 1
                if leng >= len(seq):  # have read enough quality
                    last = None
                    yield name, seq, ''.join(seqs)  # yield a fastq record
                    break
            if last:  # reach EOF before reading enough quality
                yield name, seq, None  # yield a fasta record instead
                break


def case_convert(seq):
    base_code = {'Y': ['C', 'T'], 'R': ['A', 'G'], 'W': ['A', 'T'], 'S': ['G', 'C'], 'K': ['T', 'G'], 'M': ['C', 'A'],
                 'D': ['A', 'G', 'T'], 'V': ['A', 'C', 'G'], 'H': ['A', 'C', 'T'], 'B': ['C', 'G', 'T'],
                 'N': ['A', 'T', 'C', 'G'], 'X': ['A', 'T', 'C', 'G']}

    up_string = seq.upper()
    up_list = list(up_string)
    for i in xrange(len(up_list)):
        if up_list[i] in base_code:
            up_list[i] = random.choice(base_code[up_list[i]])
    out_seq = ''.join(up_list)

    return out_seq


def simulation_aligned_transcriptome(model_ir, out_reads, out_error, kmer_bias, basecaller, read_type, num_simulate,
                                     per=False, uracil=False):
    # Simulate aligned reads
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    if model_ir:
        # check whether chrom names contains "chr" or not.
        flag_chrom = False
        for item in genome_fai.references:
            if "chr" in item:
                flag_chrom = True
                break

    sampled_2d_lengths = get_length_kde(kde_aligned_2d, num_simulate, False, False)

    remainder_l = get_length_kde(kde_ht, num_simulate, True)
    head_vs_ht_ratio_temp = get_length_kde(kde_ht_ratio, num_simulate)
    head_vs_ht_ratio_l = [1 if x > 1 else x for x in head_vs_ht_ratio_temp]
    head_vs_ht_ratio_l = [0 if x < 0 else x for x in head_vs_ht_ratio_l]

    i = 0
    while i < num_simulate:
        while True:
            ref_trx, ref_trx_len = select_ref_transcript(ecdf_dict_ref_exp)
            if model_ir:
                if ref_trx in dict_ref_structure:
                    ref_trx_len_fromstructure = ref_len_from_structure(dict_ref_structure[ref_trx])
                    if ref_trx_len == ref_trx_len_fromstructure:
                        ref_len_aligned = select_nearest_kde2d(sampled_2d_lengths, ref_trx_len)
                        if ref_len_aligned < ref_trx_len:
                            break
            else:
                ref_len_aligned = select_nearest_kde2d(sampled_2d_lengths, ref_trx_len)
                if ref_len_aligned < ref_trx_len:
                    break
        if per:
            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            new_read, ref_start_pos = extract_read_trx(ref_trx, ref_len_aligned)
            new_read_name = ref_trx + "_" + str(ref_start_pos) + "_perfect_" + str(sequence_index)
            read_mutated = case_convert(new_read)  # not mutated actually, just to be consistent with per == False
        else:
            middle_read, middle_ref, error_dict = error_list(ref_len_aligned, match_markov_model, match_ht_list,
                                                             error_par, trans_error_pr)

            if middle_ref > ref_trx_len:
                continue

            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            if model_ir:
                ir_flag, ref_trx_structure_new = update_structure(dict_ref_structure[ref_trx], IR_markov_model)
                if ir_flag:
                    list_iv = extract_read_pos(middle_ref, ref_trx_len, ref_trx_structure_new)
                    new_read = ""
                    flag = False
                    for interval in list_iv:
                        chrom = interval.chrom
                        if flag_chrom:
                            chrom = "chr" + chrom
                        if chrom not in genome_fai.references:
                            flag = True
                            break
                        start = interval.start
                        end = interval.end
                        new_read += genome_fai.fetch(chrom, start, end)
                    if flag:
                        continue
                    ref_start_pos = list_iv[0].start
                else:
                    new_read, ref_start_pos = extract_read_trx(ref_trx, middle_ref)

            else:
                new_read, ref_start_pos = extract_read_trx(ref_trx, middle_ref)

            new_read_name = str(ref_trx) + "_" + str(ref_start_pos) + "_aligned_" + str(sequence_index)

            # start HD len simulation
            remainder = int(remainder_l[i])
            head_vs_ht_ratio = head_vs_ht_ratio_l[i]

            if remainder == 0:
                head = 0
                tail = 0
            else:
                head = int(round(remainder * head_vs_ht_ratio))
                tail = remainder - head
            # end HD len simulation

            # Mutate read
            new_read = case_convert(new_read)
            read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)
            if kmer_bias:
                read_mutated = mutate_homo(read_mutated, kmer_bias, basecaller, read_type)

        # Reverse complement according to strandness rate
        p = random.random()
        if p < strandness_rate:
            read_mutated = reverse_complement(read_mutated)
            new_read_name += "_R"
        else:
            new_read_name += "_F"

        if per:
            out_reads.write(">" + new_read_name + "_0_" + str(ref_len_aligned) + "_0" + '\n')
        else:
            read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                           ''.join(np.random.choice(BASES, tail))

            out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle_ref) + "_" + str(tail) + '\n')

        if uracil:
            read_mutated = read_mutated.translate(trantab)
        out_reads.write(read_mutated + '\n')

        if (sequence_index + 1) % 100 == 0:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " +
                             str(sequence_index + 1) + "\r")
            # +1 is just to ignore the zero index by python
            sys.stdout.flush()

        i += 1

    out_reads.close()
    out_error.close()


def simulation_aligned_genome(dna_type, min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias, basecaller,
                              num_simulate, per=False):
        # Simulate aligned reads
        out_reads = open(out_reads, "w")
        out_error = open(out_error, "w")

        i = num_simulate
        passed = 0
        while i > 0:
            if per:
                ref_l = get_length_kde(kde_aligned, i) if median_l is None else \
                    np.random.lognormal(np.log(median_l), sd_l, i)
                ref_l = [x for x in ref_l if min_l <= x <= max_l and x <= max_chrom]
            else:
                remainder_l = get_length_kde(kde_ht, i, True)
                head_vs_ht_ratio_l = get_length_kde(kde_ht_ratio, i)
                if median_l is None:
                    ref_l = get_length_kde(kde_aligned, i)
                else:
                    total_l = np.random.lognormal(np.log(median_l + sd_l ** 2 / 2), sd_l, i)
                    ref_l = total_l - remainder_l

                ref_l = [x for x in ref_l if x > 0]

            for j in xrange(len(ref_l)):
                # check if the total length fits the criteria
                ref = int(ref_l[j])

                if per:
                    with total_simulated.get_lock():
                        sequence_index = total_simulated.value
                        total_simulated.value += 1

                    # Extract middle region from reference genome
                    new_read, new_read_name = extract_read(dna_type, ref)
                    new_read_name = new_read_name + "_perfect_" + str(sequence_index)
                    read_mutated = case_convert(
                        new_read)  # not mutated actually, just to be consistent with per == False
                else:
                    middle, middle_ref, error_dict = error_list(ref, match_markov_model, match_ht_list, error_par,
                                                                trans_error_pr)
                    remainder = int(remainder_l[j])
                    head_vs_ht_ratio = head_vs_ht_ratio_l[j]

                    total = remainder + middle

                    if total < min_l or total > max_l or head_vs_ht_ratio < 0 or head_vs_ht_ratio > 1 or \
                            middle_ref > max_chrom or total > max_chrom:
                        continue

                    with total_simulated.get_lock():
                        sequence_index = total_simulated.value
                        total_simulated.value += 1

                    if remainder == 0:
                        head = 0
                        tail = 0
                    else:
                        head = int(round(remainder * head_vs_ht_ratio))
                        tail = remainder - head

                    # Extract middle region from reference genome
                    new_read, new_read_name = extract_read(dna_type, middle_ref)
                    new_read_name = new_read_name + "_aligned_" + str(sequence_index)

                    # Mutate read
                    new_read = case_convert(new_read)
                    read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)
                    if kmer_bias:
                        read_mutated = mutate_homo(read_mutated, kmer_bias, basecaller, None)

                # Reverse complement half of the reads
                p = random.random()
                if p < strandness_rate:
                    read_mutated = reverse_complement(read_mutated)
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"

                if per:
                    out_reads.write(">" + new_read_name + "_0_" + str(ref) + "_0" + '\n')
                else:
                    # Add head and tail region
                    read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                                   ''.join(np.random.choice(BASES, tail))

                    out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle_ref) + "_" +
                                    str(tail) + '\n')
                out_reads.write(read_mutated + '\n')

                if (sequence_index + 1) % 100 == 0:
                    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " +
                                     str(sequence_index + 1) + "\r")
                    # +1 is just to ignore the zero index by python
                    sys.stdout.flush()

                passed += 1
            i = num_simulate - passed

        out_reads.close()
        out_error.close()


def simulation_unaligned(dna_type, min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias, basecaller, read_type,
                         num_simulate, uracil):
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    i = num_simulate
    passed = 0

    while i > 0:
        # if the median length and sd is set, use log normal distribution for simulation
        ref_l = get_length_kde(kde_unaligned, i) if median_l is None else \
            np.random.lognormal(np.log(median_l), sd_l, i)

        for j in xrange(len(ref_l)):
            # check if the total length fits the criteria
            ref = int(ref_l[j])

            unaligned, middle_ref, error_dict = unaligned_error_list(ref, error_par)
            if unaligned < min_l or unaligned > max_l or middle_ref > max_chrom:
                continue

            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            new_read, new_read_name = extract_read(dna_type, middle_ref)
            new_read_name = new_read_name + "_unaligned_" + str(sequence_index)
            # Change lowercase to uppercase and replace N with any base
            new_read = case_convert(new_read)
            read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)
            if kmer_bias:
                read_mutated = mutate_homo(read_mutated, kmer_bias, basecaller, read_type)

            # Reverse complement some of the reads based on direction information
            p = random.random()
            if p < strandness_rate:
                read_mutated = reverse_complement(read_mutated)
                new_read_name += "_R"
            else:
                new_read_name += "_F"

            out_reads.write(">" + new_read_name + "_0_" + str(middle_ref) + "_0" + '\n')
            if uracil:
                read_mutated = read_mutated.traslate(trantab)
            out_reads.write(read_mutated + "\n")

            if (sequence_index + 1) % 100 == 0:
                sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " +
                                 str(sequence_index + 1) + "\r")
                # +1 is just to ignore the zero index by python
                sys.stdout.flush()
            passed += 1

        i = num_simulate - passed

    out_reads.close()
    out_error.close()


def simulation(mode, out, dna_type, per, kmer_bias, basecaller, read_type, max_l, min_l, num_threads, median_l=None,
               sd_l=None, model_ir=False, uracil=False):
    global total_simulated  # Keeps track of number of reads that have been simulated so far
    total_simulated = mp.Value("i", 0, lock=True)

    # Start simulation
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    sys.stdout.flush()

    if num_threads == 1:
        out_aligned = out + "_aligned_reads.fasta"
        out_unaligned = out + "_unaligned_reads.fasta"
        out_error_aligned = out + "_aligned_error_profile"
        out_error_unaligned = out + "_unaligned_error_profile"
        if mode == "genome":
            simulation_aligned_genome(dna_type, min_l, max_l, median_l, sd_l, out_aligned, out_error_aligned,
                                      kmer_bias, basecaller, number_aligned, per)
        else:
            simulation_aligned_transcriptome(model_ir, out_aligned, out_error_aligned, kmer_bias, basecaller,
                                             read_type, number_aligned, per, uracil)

        if not per:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of random reads\n")
            sys.stdout.flush()

            simulation_unaligned(dna_type, min_l, max_l, median_l, sd_l, out_unaligned, out_error_unaligned,
                                 kmer_bias, basecaller, read_type, number_unaligned, uracil)

        with open(out + "_error_profile", "w") as out_error:
            with open(out_error_aligned) as infile:
                out_error.write(infile.read())

            with open(out_error_unaligned) as infile:
                out_error.write(infile.read())
        os.remove(out_error_aligned)
        os.remove(out_error_unaligned)

    else:
        procs = []
        aligned_subfiles = []
        error_subfiles = []
        num_simulate = int(number_aligned / num_threads)

        if mode == "genome":
            for i in range(num_threads):
                aligned_subfile = out + "_aligned_reads{}.fasta".format(i)
                error_subfile = out + "_error_profile{}.fasta".format(i)
                aligned_subfiles.append(aligned_subfile)
                error_subfiles.append(error_subfile)

                if i != num_threads - 1:
                    p = mp.Process(target=simulation_aligned_genome, args=(dna_type, min_l, max_l, median_l, sd_l,
                                                                           aligned_subfile, error_subfile, kmer_bias,
                                                                           basecaller, num_simulate, per))
                    procs.append(p)
                    p.start()
                else:  # Last process will simulate the remaining reads
                    p = mp.Process(target=simulation_aligned_genome, args=(dna_type, min_l, max_l, median_l, sd_l,
                                                                           aligned_subfile, error_subfile, kmer_bias,
                                                                           basecaller,
                                                                           num_simulate + number_aligned % num_threads,
                                                                           per))
                    procs.append(p)
                    p.start()

            for p in procs:
                p.join()

        else:
            for i in range(num_threads):
                aligned_subfile = out + "_aligned_reads{}.fasta".format(i)
                error_subfile = out + "_error_profile{}.fasta".format(i)
                aligned_subfiles.append(aligned_subfile)
                error_subfiles.append(error_subfile)

                if i != num_threads - 1:
                    p = mp.Process(target=simulation_aligned_transcriptome, args=(model_ir, aligned_subfile,
                                                                                  error_subfile, kmer_bias, basecaller,
                                                                                  read_type, num_simulate, per, uracil))
                    procs.append(p)
                    p.start()
                else:
                    p = mp.Process(target=simulation_aligned_transcriptome, args=(model_ir, aligned_subfile, error_subfile,
                                                                                  kmer_bias, basecaller, read_type,
                                                                                  num_simulate + number_aligned % num_threads,
                                                                                  per, uracil))
                    procs.append(p)
                    p.start()

            for p in procs:
                p.join()

        # Merging aligned reads subfiles
        with open(out + "_aligned_reads.fasta", 'w') as out_aligned_reads:
            for fname in aligned_subfiles:
                with open(fname) as infile:
                    out_aligned_reads.write(infile.read())

        for fname in aligned_subfiles:
            os.remove(fname)

        # Simulate unaligned reads, if per, number_unaligned = 0, taken care of in read_ecdf
        if not per:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of random reads\n")
            sys.stdout.flush()
            unaligned_subfiles = []
            num_simulate = int(number_unaligned / num_threads)
            for i in range(num_threads):
                unaligned_subfile = out + "_unaligned_reads{}.fasta".format(i)
                # Named "num_threads + i" so file name does not overlap with error files from aligned reads
                error_subfile = out + "_error_profile{}.fasta".format(num_threads + i)
                unaligned_subfiles.append(unaligned_subfile)
                error_subfiles.append(error_subfile)

                # Dividing number of unaligned reads that need to be simulated amongst the number of processes
                if i != num_threads - 1:
                    p = mp.Process(target=simulation_unaligned, args=(dna_type, min_l, max_l, median_l, sd_l,
                                                                      unaligned_subfile, error_subfile, kmer_bias,
                                                                      basecaller, read_type, num_simulate, uracil))
                    procs.append(p)
                    p.start()
                else:
                    p = mp.Process(target=simulation_unaligned, args=(dna_type, min_l, max_l, median_l, sd_l,
                                                                      unaligned_subfile, error_subfile, kmer_bias,
                                                                      basecaller, read_type,
                                                                      num_simulate + number_unaligned % num_threads,
                                                                      uracil))
                    procs.append(p)
                    p.start()

            for p in procs:
                p.join()

            # Merging unaligned reads subfiles
            with open(out + "_unaligned_reads.fasta", 'w') as out_unaligned_reads:
                for fname in unaligned_subfiles:
                    with open(fname) as infile:
                        out_unaligned_reads.write(infile.read())

            for fname in unaligned_subfiles:
                os.remove(fname)

        # Merging error subfiles
        with open(out + "_error_profile", 'w') as out_error:
            out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")
            for fname in error_subfiles:
                with open(fname) as infile:
                    out_error.write(infile.read())

        for fname in error_subfiles:
            os.remove(fname)


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def extract_read_trx(key, length):
    ref_pos = random.randint(0, seq_len[key] - length)
    new_read = seq_dict[key][ref_pos: ref_pos + length]
    return new_read, ref_pos


def extract_read(dna_type, length):
    if dna_type == "transcriptome":
        while True:
            new_read = ""
            key = random.choice(list(seq_len.keys()))  # added "list" thing to be compatible with Python v3
            if length < seq_len[key]:
                ref_pos = random.randint(0, seq_len[key] - length)
                new_read = seq_dict[key][ref_pos: ref_pos + length]
                new_read_name = key + "_" + str(ref_pos)
                break
        return new_read, new_read_name
    else:
        # Extract the aligned region from reference
        if dna_type == "circular":
            ref_pos = random.randint(0, genome_len)
            chromosome = list(seq_dict.keys())[0]
            new_read_name = chromosome + "_" + str(ref_pos)
            if length + ref_pos <= genome_len:
                new_read = seq_dict[chromosome][ref_pos: ref_pos + length]
            else:
                new_read = seq_dict[chromosome][ref_pos:]
                new_read = new_read + seq_dict[chromosome][0: length - genome_len + ref_pos]
        else:
            # Generate a random number within the size of the genome. Suppose chromosomes are connected
            # tail to head one by one in the order of the dictionary. If the start position fits in one
            # chromosome, but the end position does not, then restart generating random number.
            while True:
                new_read = ""
                ref_pos = random.randint(0, genome_len)
                for key in seq_len.keys():
                    if ref_pos + length <= seq_len[key]:
                        new_read = seq_dict[key][ref_pos: ref_pos + length]
                        new_read_name = key + "_" + str(ref_pos)
                        break
                    elif ref_pos < seq_len[key]:
                        break
                    else:
                        ref_pos -= seq_len[key]
                if new_read != "":
                    break
        return new_read, new_read_name


def unaligned_error_list(m_ref, error_p):
    l_new = m_ref
    e_dict = {}
    error_rate = {(0, 0.4): "match", (0.4, 0.7): "mis", (0.7, 0.85): "ins", (0.85, 1): "del"}
    pos = 0
    middle_ref = m_ref
    last_is_ins = False
    while pos < middle_ref:
        p = random.random()
        for k_error in error_rate.keys():
            if k_error[0] <= p < k_error[1]:
                error_type = error_rate[k_error]
                break

        if error_type == "match":
            step = 1

        elif error_type == "mis":
            step = mm.pois_geom(error_p["mis"][0], error_p["mis"][2], error_p["mis"][3])
            e_dict[pos] = ["mis", step]

        elif error_type == "ins":
            step = mm.wei_geom(error_p["ins"][0], error_p["ins"][1], error_p["ins"][2], error_p["ins"][3])
            if last_is_ins:
                e_dict[pos + 0.1][1] += step
            else:
                e_dict[pos + 0.1] = ["ins", step]
                last_is_ins = True
            l_new += step

        else:
            step = mm.wei_geom(error_p["del"][0], error_p["del"][1], error_p["del"][2], error_p["del"][3])
            e_dict[pos] = ["del", step]
            l_new -= step

        if error_type != "ins":
            pos += step
            last_is_ins = False

        if pos > middle_ref:
            l_new += pos - middle_ref
            middle_ref = pos

    return l_new, middle_ref, e_dict


def error_list(m_ref, m_model, m_ht_list, error_p, trans_p):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    middle_ref = m_ref
    prev_error = "start"

    # The first match come from m_ht_list
    p = random.random()
    k1 = list(m_ht_list.keys())[0]
    for k2, v2 in m_ht_list[k1].items():
        if k2[0] < p <= k2[1]:
            prev_match = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
            if prev_match < 2:
                prev_match = 2
    pos += prev_match

    # Select an error, then the step size, and then a match and so on so forth.
    while pos < middle_ref:
        # pick the error based on Markov chain
        p = random.random()
        for k in trans_p[prev_error].keys():
            if k[0] <= p < k[1]:
                error = trans_p[prev_error][k]
                break

        if error == "mis":
            step = mm.pois_geom(error_p[error][0], error_p[error][2], error_p[error][3])
        elif error == "ins":
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new += step
        else:
            step = mm.wei_geom(error_p[error][0], error_p[error][1], error_p[error][2], error_p[error][3])
            l_new -= step

        if error != "ins":
            e_dict[pos] = [error, step]
            pos += step
            if pos >= middle_ref:
                l_new += pos - middle_ref
                middle_ref = pos
        else:
            e_dict[pos - 0.5] = [error, step]

        prev_error = error

        # Randomly select a match length
        for k1 in m_model.keys():
            if k1[0] <= prev_match < k1[1]:
                break
        p = random.random()
        for k2, v2 in m_model[k1].items():
            if k2[0] < p <= k2[1]:
                step = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
                break
        # there are no two 0 base matches together
        if prev_match == 0 and step == 0:
            step = 1

        prev_match = step
        if pos + prev_match > middle_ref:
            l_new += pos + prev_match - middle_ref
            middle_ref = pos + prev_match

        pos += prev_match
        if prev_match == 0:
            prev_error += "0"

    return l_new, middle_ref, e_dict


def mutate_read(read, read_name, error_log, e_dict, k, aligned=True):
    if k:  # First remove any errors that land in hp regions
        pattern = "A{" + re.escape(str(k)) + ",}|C{" + re.escape(str(k)) + ",}|G{" + re.escape(str(k)) + ",}|T{" + \
                  re.escape(str(k)) + ",}"

        hp_pos = []  # [[start, end], ...]
        for match in re.finditer(pattern, read):
            hp_pos.append([match.start(), match.end()])

        new_e_dict = {}
        for err_start in e_dict.keys():
            err = e_dict[err_start][0]
            err_end = err_start + e_dict[err_start][1]
            hp_err = False

            for hp in hp_pos:
                hp_start = hp[0]
                hp_end = hp[1]
                if not (hp_end <= err_start or err_end <= hp_start):  # Lands in hp
                    hp_err = True
                    break

            if not hp_err:
                new_e_dict[err_start] = [err, e_dict[err_start][1]]
    else:
        new_e_dict = e_dict

    # Mutate read
    for key in sorted(new_e_dict.keys(), reverse=True):
        val = new_e_dict[key]
        key = int(round(key))

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            new_bases = ""
            for i in xrange(val[1]):
                tmp_bases = list(BASES)
                tmp_bases.remove(read[key + i])
                # tmp_bases.remove(read[key]) ## Edited this part for testing
                new_base = random.choice(tmp_bases)
                new_bases += new_base

            new_read = read[:key] + new_bases + read[key + val[1]:]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            new_bases = ""
            for i in xrange(val[1]):
                new_base = random.choice(BASES)
                new_bases += new_base
            new_read = read[:key] + new_bases + read[key:]

        read = new_read

        if aligned and val[0] != "match":
            error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            "\t" + ref_base + "\t" + new_bases + "\n")

    return read


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Simulation step
        -----------------------------------------------------------
        Given error profiles, reference genome and/or transcriptome,
        simulate ONT DNA or RNA reads
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='NanoSim ' + VERSION)
    subparsers = parser.add_subparsers(help="You may run the simulator on transcriptome or genome mode.", dest='mode',
                                       description=dedent('''
        There are two modes in read_analysis.
        For detailed usage of each mode:
            simulator.py mode -h
        -------------------------------------------------------
        '''))

    parser_g = subparsers.add_parser('genome', help="Run the simulator on genome mode")
    parser_g.add_argument('-rg', '--ref_g', help='Input reference genome', required=True)
    parser_g.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                       'characterization step (Default = training)',
                          default="training")
    parser_g.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                          default="simulated")
    parser_g.add_argument('-n', '--number', help='Number of reads to be simulated (Default = 20000)', type=int,
                          default=20000)
    parser_g.add_argument('-max', '--max_len', help='The maximum length for simulated reads (Default = Infinity)',
                          type=int, default=float("inf"))
    parser_g.add_argument('-min', '--min_len', help='The minimum length for simulated reads (Default = 50)',
                          type=int, default=50)
    parser_g.add_argument('-med', '--median_len', help='The median read length (Default = None)',
                          type=int, default=None)
    parser_g.add_argument('-sd', '--sd_len', help='The standard deviation of read length in log scale (Default = None)',
                          type=float, default=None)
    parser_g.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_g.add_argument('-k', '--KmerBias', help='Minimum homopolymer length to simulate homopolymer contraction and '
                                                   'expansion events in',
                          type=int, default=None)
    parser_g.add_argument('-b', '--basecaller', help='Simulate homopolymers with respect to chosen basecaller: '
                                                     'albacore, guppy, or guppy-flipflop',
                          choices=["albacore", "guppy", "guppy-flipflop"], default=None)
    parser_g.add_argument('-s', '--strandness', help='Percentage of antisense sequences. Overrides the value profiled '
                                                     'in characterization stage. Should be between 0 and 1',
                          type=float, default=None)
    parser_g.add_argument('-dna_type', help='Specify the dna type: circular OR linear (Default = linear)',
                          choices=["linear", "circular"], default="linear")
    parser_g.add_argument('--perfect', help='Ignore error profiles and simulate perfect reads', action='store_true',
                          default=False)
    parser_g.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int, default=1)

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode")
    parser_t.add_argument('-rt', '--ref_t', help='Input reference transcriptome', required=True)
    parser_t.add_argument('-rg', '--ref_g', help='Input reference genome, required if intron retention simulatin is on',
                          default='')
    parser_t.add_argument('-e', '--exp', help='Expression profile in the specified format as described in README',
                          required=True)
    parser_t.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                       'characterization step (Default = training)',
                          default="training")
    parser_t.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                          default="simulated")
    parser_t.add_argument('-n', '--number', help='Number of reads to be simulated (Default = 20000)', type=int,
                          default=20000)
    parser_t.add_argument('-max', '--max_len', help='The maximum length for simulated reads (Default = Infinity)',
                          type=int, default=float("inf"))
    parser_t.add_argument('-min', '--min_len', help='The minimum length for simulated reads (Default = 50)',
                          type=int, default=50)
    parser_t.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_t.add_argument('-k', '--KmerBias', help='Minimum homopolymer length to simulate homopolymer contraction and '
                                                   'expansion events in',
                          type=int, default=None)
    parser_t.add_argument('-b', '--basecaller', help='Simulate homopolymers with respect to chosen basecaller: '
                                                     'albacore or guppy',
                          choices=["albacore", "guppy"], default=None)
    parser_t.add_argument('-r', '--read_type', help='Simulate homopolymers with respect to chosen read type: dRNA, '
                                                    'cDNA_1D or cDNA_1D2',
                          choices=["dRNA", "cDNA_1D", "cDNA_1D2"], default=None)
    parser_t.add_argument('-s', '--strandness', help='Percentage of antisense sequences. Overrides the value profiled '
                                                     'in characterization stage. Should be between 0 and 1',
                          type=float, default=None)
    parser_t.add_argument('--no_model_ir', help='Ignore simulating intron retention events', action='store_false', default=True)
    parser_t.add_argument('--perfect', help='Ignore profiles and simulate perfect reads', action='store_true',
                          default=False)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int, default=1)
    parser_t.add_argument('--uracil', help='Converts the thymine (T) bases to uracil (U) in the output fasta format', action='store_true', default=False)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.mode == "genome":
        ref_g = args.ref_g
        model_prefix = args.model_prefix
        out = args.output
        number = args.number
        max_len = args.max_len
        min_len = args.min_len
        median_len = args.median_len
        sd_len = args.sd_len
        if args.seed:
            random.seed(int(args.seed))
            np.random.seed(int(args.seed))
        perfect = args.perfect
        kmer_bias = args.KmerBias
        basecaller = args.basecaller
        strandness = args.strandness
        dna_type = args.dna_type
        num_threads = max(args.num_threads, 1)

        if kmer_bias and kmer_bias < 0:
            print("\nPlease input proper kmer bias value >= 0\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if kmer_bias and basecaller is None:
            print("\nPlease input basecaller to simulate homopolymer contraction and expansion events from\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if strandness and (strandness < 0 or strandness > 1):
            print("\nPlease input proper strandness value between 0 and 1\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if (median_len and not sd_len) or (sd_len and not median_len):
            sys.stderr.write("\nPlease provide both mean and standard deviation of read length!\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if max_len < min_len:
            sys.stderr.write("\nMaximum read length must be longer than Minimum read length!\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("ref_g", ref_g)
        print("model_prefix", model_prefix)
        print("out", out)
        print("number", number)
        print("perfect", perfect)
        print("kmer_bias", kmer_bias)
        print("basecaller", basecaller)
        print("dna_type", dna_type)
        print("strandness", strandness)
        print("sd_len", sd_len)
        print("median_len", median_len)
        print("max_len", max_len)
        print("min_len", min_len)
        print("num_threads", num_threads)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        read_profile(ref_g, None, number, model_prefix, perfect, args.mode, strandness, None, False, dna_type)

        if median_len and sd_len:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Simulating read length with log-normal distribution\n")
            sys.stdout.flush()
            simulation(args.mode, out, dna_type, perfect, kmer_bias, basecaller, None, max_len, min_len, num_threads,
                       median_len, sd_len)
        else:
            simulation(args.mode, out, dna_type, perfect, kmer_bias, basecaller, None, max_len, min_len, num_threads)

    elif args.mode == "transcriptome":
        ref_g = args.ref_g
        ref_t = args.ref_t
        exp = args.exp
        model_prefix = args.model_prefix
        out = args.output
        number = args.number
        max_len = args.max_len
        min_len = args.min_len
        kmer_bias = args.KmerBias
        basecaller = args.basecaller
        read_type = args.read_type
        strandness = args.strandness
        perfect = args.perfect
        model_ir = args.no_model_ir
        dna_type = "transcriptome"
        uracil = args.uracil
        num_threads = max(args.num_threads, 1)

        if kmer_bias and kmer_bias < 0:
            print("\nPlease input proper kmer bias value >= 0\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        if kmer_bias and (basecaller is None or read_type is None):
            print("\nPlease input basecaller and read_type to simulate homopolymer contraction and expansion events "
                  "from\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        if strandness and (strandness < 0 or strandness > 1):
            print("\nPlease input proper strandness value between 0 and 1\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        if max_len < min_len:
            sys.stderr.write("\nMaximum read length must be longer than Minimum read length!\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        if model_ir and ref_g == '':
            sys.stderr.write("\nPlease provide a reference genome to simulate intron retention events!\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("ref_g", ref_g)
        print("ref_t", ref_t)
        print("exp", exp)
        print("model_prefix", model_prefix)
        print("out", out)
        print("number", number)
        print("perfect", perfect)
        print("kmer_bias", kmer_bias)
        print("basecaller", basecaller)
        print("read_type", read_type)
        print("model_ir", model_ir)
        print("dna_type", dna_type)
        print("strandness", strandness)
        print("max_len", max_len)
        print("min_len", min_len)
        print("uracil", uracil)
        print("num_threads", num_threads)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        read_profile(ref_g, ref_t, number, model_prefix, perfect, args.mode, strandness, exp, model_ir, "linear")

        simulation(args.mode, out, dna_type, perfect, kmer_bias, basecaller, read_type, max_len, min_len, num_threads,
                   None, None, model_ir, uracil)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
    sys.stdout.close()


if __name__ == "__main__":
    main()
