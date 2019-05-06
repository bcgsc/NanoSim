#!/usr/bin/env python
"""

@author: Chen Yang & Saber HafezQorani

This script generates simulated Oxford Nanopore 2D reads (genomic and transcriptomic - cDNA/directRNA).

"""


from __future__ import print_function
from __future__ import with_statement
from subprocess import call
import sys
import os
import HTSeq
import pysam
import numpy
import random
import re
import copy
import argparse
from time import strftime
from time import sleep
import numpy as np
from sklearn.externals import joblib
#from math import exp

try:
    from six.moves import xrange
except ImportError:
    pass
import mixed_model as mm

PYTHON_VERSION = sys.version_info
VERSION = "2.1.0"
PRORAM = "NanoSim"
AUTHOR = "Chen Yang (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca"

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
    for i in range (0, len(l) - 1):
        r = (l[i], l[i+1])
        output_list.append(r)
    return output_list


def make_cdf(dict_exp, dict_len):
    sum_exp = 0
    list_value = []
    dict_tpm = {}
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
        #ecdf_dict[cdf_range] = dict_len[sorted_value_list[i][0]]

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
    prev_state = "start"
    for i in range (0, count):
        p = random.random()
        for key in IR_markov_model[prev_state]:
            if key[0] <= p < key[1]:
                flag = IR_markov_model[prev_state][key]
                list_states.append(flag)
                prev_state = flag

    j = -1
    for i in range (0, len(ref_trx_structure)):
        if ref_trx_structure[i][0] == "intron":
            j += 1
            if list_states[j] == "IR":
                ref_trx_structure[i] = ("retained_intron",) + ref_trx_structure[i][1:]

    return list_states, ref_trx_structure


def extract_read_pos(length, ref_len, ref_trx_structure):
    #The aim is to create a genomic interval object
    #example: iv = HTSeq.GenomicInterval( "chr3", 123203, 127245, "+" )
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
            if flag == False: #if it is first exon that I started extracting from
                if start_pos < item[-1]:
                    #print ("1", item, start_pos, length, start, end)
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
                    #print("2", item, start_pos, length, start, end)
                    start_pos -= item[-1]
            else: #if it is NOT the first exon that I start extracting from
                #print("3", item, start_pos, length, start, end)
                start_pos = 0
                if (item[2] + length) < item[3]:
                    end = item[2] + length
                    length = 0

                else:
                    end = item[3]
                    length -= (item[-1])

        elif item[0] == "retained_intron":
            #print("4", item, start_pos, length, start, end)
            if flag != False:
                end = item[3]
        elif item[0] == "intron":
            #print("5", item, start_pos, length, start, end)
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


def read_profile(ref_g, ref_t, number, model_prefix, per, mode, strandness, exp=None, model_ir=None):
    global number_aligned, number_unaligned
    global match_ht_list, error_par, trans_error_pr, match_markov_model
    global kde_aligned, kde_ht, kde_ht_ratio, kde_unaligned, kde_aligned_2d
    global seq_dict, seq_len
    global strandness_rate

    if mode == "genome":
        global genome_len
        ref = ref_g
    else:
        global dict_ref_structure, dict_exp, ecdf_dict_ref_exp, genome_fai
        ref = ref_t

    if strandness == None:
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

    if mode == "genome":
        genome_len = sum(seq_len.values())
        if len(seq_dict) > 1 and dna_type == "circular":
            sys.stderr.write("Do not choose circular if there is more than one chromosome in the genome!")
            sys.exit(1)
    else:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference genome and create .fai index file\n")
        sys.stdout.flush()
        # create and read the .fai file of the reference genome
        genome_fai = pysam.Fastafile(ref_g)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in expression profile\n")
        sys.stdout.flush()
        dict_exp = {}
        with open (exp, 'r') as exp_file:
            header = exp_file.readline()
            for line in exp_file:
                parts = line.split("\t")
                transcript_id = parts[0].split(".")[0]
                tpm = float(parts[2])
                if transcript_id.startswith("ENS") and tpm > 0:
                    dict_exp[transcript_id] = tpm

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in GFF3 annotation file\n")
        sys.stdout.flush()
        dict_ref_structure = {}
        gff_file = model_prefix + "_addedintron_final.gff3"
        gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)
        for feature in gff_features:
            if feature.type == "exon" or feature.type == "intron":
                flag = False
                if "transcript_id" in feature.attr:
                    feature_id = feature.attr['transcript_id']
                    flag = True
                elif "Parent" in feature.attr:
                    info = feature.name.split(":")
                    if len(info) == 1:
                        feature_id = info[0]
                    else:
                        if info[0] == "transcript":
                            feature_id = info[1]
                        else:
                            continue
                    flag = True

                if flag and "ENS" in feature_id:
                    feature_id = feature_id.split(".")[0]
                    if feature_id not in dict_ref_structure:
                        dict_ref_structure[feature_id] = []

                    # remove "chr" from chromosome names to be constant
                    if "chr" in feature.iv.chrom:
                        feature.iv.chrom = feature.iv.chrom.strip("chr")

                    dict_ref_structure[feature_id].append(
                        (feature.type, feature.iv.chrom, feature.iv.start, feature.iv.end, feature.iv.length))

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

    kde_unaligned = joblib.load(model_prefix + "_unaligned_length.pkl")

    with open(model_prefix + "_reads_alignment_rate", 'r') as u_profile:
        new = u_profile.readline().strip()
        rate = new.split('\t')[1]
        # if parameter perfect is used, all reads should be aligned, number_aligned equals total number of reads.
        if per or rate == "100%":
            number_aligned = number
        else:
            number_aligned = int(round(number * float(rate) / (float(rate) + 1)))
        number_unaligned = number - number_aligned

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


def collapse_homo(seq, k):
    read = re.sub("A" * k + "+", "A" * (k - 1), seq)
    read = re.sub("C" * k + "+", "C" * (k - 1), read)
    read = re.sub("T" * k + "+", "T" * (k - 1), read)
    read = re.sub("G" * k + "+", "G" * (k - 1), read)

    return read


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


def simulation_aligned_transcriptome(model_ir, out_reads, out_error, kmer_bias):

    # Simulate aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    sys.stdout.flush()
    i = 0
    sampled_2d_lengths = get_length_kde(kde_aligned_2d, number_aligned, False, False)
    remainder_l = get_length_kde(kde_ht, number_aligned, True)
    head_vs_ht_ratio_l = get_length_kde(kde_ht_ratio, number_aligned)
    head_vs_ht_ratio_l_new = []
    for x in head_vs_ht_ratio_l:
        if 0 < x < 1:
            head_vs_ht_ratio_l_new.append(x)
        elif x < 0:
            head_vs_ht_ratio_l_new.append(0)
        elif x > 1:
            head_vs_ht_ratio_l_new.append(1)

    while i < number_aligned:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(i + 1 + number_unaligned) + "\r") #+1 is just to ignore the zero index by python
        sys.stdout.flush()
        sleep(0.02)
        while True:
            ref_trx, ref_trx_len = select_ref_transcript(ecdf_dict_ref_exp)
            if ref_trx in dict_ref_structure:
                ref_trx_structure = copy.deepcopy(dict_ref_structure[ref_trx])
                ref_trx_len_fromstructure = ref_len_from_structure(ref_trx_structure)
                if ref_trx_len == ref_trx_len_fromstructure:
                    ref_len_aligned = select_nearest_kde2d(sampled_2d_lengths, ref_trx_len)
                    if ref_len_aligned < ref_trx_len:
                        break

        #ir_length = 0
        if model_ir:
            ir_info, ref_trx_structure_new = update_structure(ref_trx_structure, IR_markov_model)
            # for item in ref_trx_structure_new:
            #     if item[0] == "retained_intron":
            #         ir_length += item[-1]
        else:
            ref_trx_structure_new = copy.deepcopy(dict_ref_structure[ref_trx])

        #check whether chrom names contains "chr" or not.
        flag_chrom = False
        for item in genome_fai.references:
            if "chr" in item:
                flag_chrom = True

        list_iv = extract_read_pos(ref_len_aligned, ref_trx_len, ref_trx_structure_new)
        new_read = ""
        flag = False
        for interval in list_iv:
            chrom = interval.chrom
            if flag_chrom == True:
                chrom = "chr" + chrom
            if chrom not in genome_fai.references:
                flag = True
                break
            start = interval.start
            end = interval.end
            new_read += genome_fai.fetch(chrom, start, end)
        if flag == True:
            continue

        new_read_length = len(new_read)
        middle_read, middle_ref, error_dict = error_list(new_read_length, match_markov_model, match_ht_list, error_par,
                                                        trans_error_pr)
        #start HD len simulation
        remainder = int(remainder_l[i])
        head_vs_ht_ratio = head_vs_ht_ratio_l_new[i]

        total = remainder + middle_read

        if head_vs_ht_ratio < 0 or head_vs_ht_ratio > 1:
            continue

        if remainder == 0:
            head = 0
            tail = 0
        else:
            head = int(round(remainder * head_vs_ht_ratio))
            tail = remainder - head
        #end HD len simulation

        ref_start_pos = list_iv[0].start
        new_read_name = str(ref_trx) + "_" + str(ref_start_pos) + "_aligned_" + str(i + number_unaligned)
        # Mutate read
        new_read = case_convert(new_read)
        try: #because I first extract read and then introduce errors. they are several cases in which the error introduced length is larger than read length. [error in last mis error length] #dev
            read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)
        except:
            #print (new_read_length, middle_read, middle_ref, sorted(error_dict.keys(), reverse=True)[0], error_dict[sorted(error_dict.keys(), reverse=True)[0]])
            continue

        # Reverse complement accoding to strandness rate
        p = random.random()
        if p < strandness_rate:
            read_mutated = reverse_complement(read_mutated)
            new_read_name += "_R"
        else:
            new_read_name += "_F"

        # Add head and tail region
        read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                       ''.join(np.random.choice(BASES, tail))

        if kmer_bias:
            read_mutated = collapse_homo(read_mutated, kmer_bias)

        out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle_ref) + "_" + str(tail) + '\n')
        out_reads.write(read_mutated + '\n')

        i += 1


def simulation_aligned_genome(dna_type, min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias):

    # Simulate aligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    sys.stdout.flush()
    i = number_aligned
    passed = 0
    while i > 0:
        remainder_l = get_length_kde(kde_ht, i, True)
        head_vs_ht_ratio_l = get_length_kde(kde_ht_ratio, i)
        if median_l is None:
            ref_l = get_length_kde(kde_aligned, i)
        else:
            total_l = np.random.lognormal(np.log(median_l), sd_l, i)
            ref_l = total_l - remainder_l

        ref_l = [x for x in ref_l if x > 0]

        for j in xrange(len(ref_l)):
            # check if the total length fits the criteria
            ref = int(ref_l[j])
            middle, middle_ref, error_dict = error_list(ref, match_markov_model, match_ht_list, error_par,
                                                        trans_error_pr)
            remainder = int(remainder_l[j])
            head_vs_ht_ratio = head_vs_ht_ratio_l[j]

            total = remainder + middle

            if total < min_l or total > max_l or head_vs_ht_ratio < 0 or head_vs_ht_ratio > 1:
                continue

            if remainder == 0:
                head = 0
                tail = 0
            else:
                head = int(round(remainder * head_vs_ht_ratio))
                tail = remainder - head

            # Extract middle region from reference genome
            new_read, new_read_name = extract_read(dna_type, middle_ref)
            new_read_name = new_read_name + "_aligned_" + str(passed + number_unaligned)

            # Mutate read
            new_read = case_convert(new_read)
            read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias)

            # Reverse complement half of the reads
            p = random.random()
            if p < strandness_rate:
                read_mutated = reverse_complement(read_mutated)
                new_read_name += "_R"
            else:
                new_read_name += "_F"

            # Add head and tail region
            read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                           ''.join(np.random.choice(BASES, tail))

            if kmer_bias:
                read_mutated = collapse_homo(read_mutated, kmer_bias)

            out_reads.write(">" + new_read_name + "_" + str(head) + "_" + str(middle_ref) + "_" +
                            str(tail) + '\n')
            out_reads.write(read_mutated + '\n')

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(
                passed + 1 + number_unaligned) + "\r")  # +1 is just to ignore the zero index by python
            sys.stdout.flush()
            sleep(0.02)

            passed += 1
        i = number_aligned - passed


def simulation(mode, out, dna_type, per, kmer_bias, max_l, min_l, median_l=None, sd_l=None, model_ir=None):
    global number_aligned, number_unaligned
    global match_ht_list, error_par, trans_error_pr, match_markov_model
    global kde_aligned, kde_ht, kde_ht_ratio, kde_unaligned, kde_aligned_2d
    global seq_dict, seq_len
    if mode == "genome":
        global genome_len
    else:
        global dict_ref_structure, IR_markov_model, dict_exp, ecdf_dict_ref_exp

    # Start simulation
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of random reads\n")
    sys.stdout.flush()
    out_reads = open(out + "_reads.fasta", 'w')
    out_error = open(out + "_error_profile", 'w')
    out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")

    # Simulate unaligned reads
    i = number_unaligned
    unaligned_length = []
    while i > 0:
        # if the median length and sd is set, use log normal distribution for simulation
        unaligned_tmp = get_length_kde(kde_unaligned, i) if median_l is None else \
            np.random.lognormal(np.log(median_l), sd_l, i)
        unaligned_length.extend([x for x in unaligned_tmp if min_l <= x <= max_l])
        i = number_unaligned - len(unaligned_length)

    for i in xrange(number_unaligned):
        unaligned = int(unaligned_length[i])
        unaligned, error_dict = unaligned_error_list(unaligned, error_par)
        new_read, new_read_name = extract_read(dna_type, unaligned)
        new_read_name = new_read_name + "_unaligned_" + str(i)
        # Change lowercase to uppercase and replace N with any base
        new_read = case_convert(new_read)
        read_mutated = mutate_read(new_read, new_read_name, out_error, error_dict, kmer_bias, False)

        # Reverse complement half of the reads
        p = random.random()
        if p < strandness_rate:
            read_mutated = reverse_complement(read_mutated)
            new_read_name += "_R"
        else:
            new_read_name += "_F"
        out_reads.write(">" + new_read_name + "_0_" + str(unaligned) + "_0" + '\n')
        out_reads.write(read_mutated + "\n")
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(i + 1) + "\r") #+1 is just to ignore the zero index by python
        sys.stdout.flush()
        sleep(0.02)

    del unaligned_length


    if per:
        i = number_aligned
        ref_length = []
        while i > 0:
            ref_tmp = get_length_kde(kde_aligned, i) if median_l is None else \
                np.random.lognormal(np.log(median_l), sd_l, i)
            ref_length.extend([x for x in ref_tmp if min_l <= x <= max_l])
            i = number_aligned - len(ref_length)

        for i in xrange(number_aligned):
            read_len = int(ref_length[i])
            new_read, new_read_name = extract_read(dna_type, read_len)
            new_read_name = new_read_name + "_perfect_" + str(i)

            # Reverse complement half of the reads
            p = random.random()
            if p < strandness_rate:
                new_read = reverse_complement(new_read)
                new_read_name += "_R"
            else:
                new_read_name += "_F"
            out_reads.write(">" + new_read_name + "_0_" + str(read_len) + "_0" + '\n')

            # Change lowercase to uppercase and replace N with any base
            new_read = case_convert(new_read)
            out_reads.write(new_read + "\n")
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " + str(i + 1 + number_unaligned) + "\r") #+1 is just to ignore the zero index by python
            sys.stdout.flush()
            sleep(0.02)

    if mode == "genome":
        simulation_aligned_genome(dna_type, min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias)
    else:
        simulation_aligned_transcriptome(model_ir, out_reads, out_error, kmer_bias)

    out_reads.close()
    out_error.close()


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def extract_read(dna_type, length):

    if dna_type == "transcriptome":
        while True:
            new_read = ""
            key = random.choice(seq_len.keys())
            if length < seq_len[key]:
                ref_pos = random.randint(0, seq_len[key] - length)
                new_read = seq_dict[key][ref_pos: ref_pos + length]
                new_read_name = key + "_" + str(ref_pos)
                break
        return new_read, new_read_name
    else:
        if length > max(seq_len.values()):
            length = max(seq_len.values())

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


def unaligned_error_list(length, error_p):
    e_dict = {}
    error_rate = {(0, 0.4): "match", (0.4, 0.7): "mis", (0.7, 0.85): "ins", (0.85, 1): "del"}
    pos = 0
    last_is_ins = False
    while pos < length:
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

        else:
            step = mm.wei_geom(error_p["del"][0], error_p["del"][1], error_p["del"][2], error_p["del"][3])
            e_dict[pos] = ["del", step]

        if error_type != "ins":
            pos += step
            last_is_ins = False

        if pos > length:
            length = pos

    return length, e_dict


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
            prev_match = int(np.floor((p - k2[0])/(k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
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
                step = int(np.floor((p - k2[0])/(k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
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
    search_pattern = "A" * k + "+|" + "T" * k + "+|" + "C" * k + "+|" + "G" * k
    for key in sorted(e_dict.keys(), reverse=True):
        val = e_dict[key]
        key = int(round(key))

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    tmp_bases = list(BASES)
                    tmp_bases.remove(read[key + i]) ##
                    #tmp_bases.remove(read[key]) ## Edited this part for testing
                    new_base = random.choice(tmp_bases)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key + val[1]: key + val[1] + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key + val[1]:]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            while True:
                new_bases = ""
                for i in xrange(val[1]):
                    new_base = random.choice(BASES)
                    new_bases += new_base
                check_kmer = read[max(key - k + 1, 0): key] + new_bases + read[key: key + k - 1]
                if not k or not re.search(search_pattern, check_kmer):
                    break
            new_read = read[:key] + new_bases + read[key:]

        read = new_read

        #if aligned and val[0] != "match":
            #error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            #"\t" + ref_base + "\t" + new_bases + "\n")

    # If choose to have kmer bias, then need to compress homopolymers to 5-mer
    if k:
        read = collapse_homo(read, k)

    return read


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


def main():

    exp = ""
    ref = ""
    model_prefix = "training"
    out = "simulated"
    number = 20000
    perfect = False
    model_ir = True
    # ins, del, mis rate represent the weight tuning in mix model
    ins_rate = 1
    del_rate = 1
    mis_rate = 1
    max_readlength = float("inf")
    min_readlength = 50
    median_readlength = None
    sd_readlength = None
    kmer_bias = 0
    strandness = None

    parser = argparse.ArgumentParser(
        description='Given the read profiles from characterization step, ' \
                    'simulate genomeic/transcriptomic ONT reads and outputs the error profiles.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    subparsers = parser.add_subparsers(help = "You may run the simulator on transcriptome or genome mode.", dest='mode')

    parser_g = subparsers.add_parser('genome', help="Run the simulator on genome mode.")
    parser_g.add_argument('-rg', '--ref_g', help='Input reference genome', type=str, required=True)
    parser_g.add_argument('-c', '--model_prefix', help='Address for profiles created in characterization step (model_prefix)', type = str, default= "training")
    parser_g.add_argument('-o', '--output', help='Output address for simulated reads', type = str, default= "simulated")
    parser_g.add_argument('-n', '--number', help='Number of reads to be simulated', type = int, default = 20000)
    parser_g.add_argument('-i', '--insertion_rate', help='Insertion rate (optional)', type = float, default= 1)
    parser_g.add_argument('-d', '--deletion_rate', help='Deletion rate (optional)', type = float, default= 1)
    parser_g.add_argument('-m', '--mismatch_rate', help='Mismatch rate (optional)', type = float, default= 1)
    parser_g.add_argument('-max', '--max_len', help='The maximum length for simulated reads', type=int, default= float("inf"))
    parser_g.add_argument('-min', '--min_len', help='The minimum length for simulated reads', type=int, default= 50)
    parser_g.add_argument('-med', '--median_len', help='The median read length', type=int, default=None)
    parser_g.add_argument('-sd', '--sd_len', help='The standard deviation of read length in log scale', type=float, default=None)
    parser_g.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_g.add_argument('-k', '--KmerBias', help='Determine whether to considert Kmer Bias or not', type = int, default= 0)
    parser_g.add_argument('-s', '--strandness', help='Determine the strandness of the simulated reads. Overrides the value profiled in characterization phase. Should be between 0 and 1.', type=float, default=None)
    parser_g.add_argument('-dna_type', help='Specify the dna type: circular OR linear, default = linear', type=str, default="linear")
    parser_g.add_argument('--perfect', help='Ignore profiles and simulate perfect reads', action='store_true')

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode.")
    parser_t.add_argument('-rt', '--ref_t', help='Input reference transcriptome', type = str, required= True)
    parser_t.add_argument('-rg', '--ref_g', help='Input reference genome', type=str, required=True)
    parser_t.add_argument('-e', '--exp', help='Expression profile in the specified format specified in the documentation', type = str)
    parser_t.add_argument('-c', '--model_prefix', help='Address for profiles created in characterization step (model_prefix)', type = str, default= "training")
    parser_t.add_argument('-o', '--output', help='Output address for simulated reads', type = str, default= "simulated")
    parser_t.add_argument('-n', '--number', help='Number of reads to be simulated', type = int, default = 20000)
    parser_t.add_argument('-i', '--insertion_rate', help='Insertion rate (optional)', type = float, default= 1)
    parser_t.add_argument('-d', '--deletion_rate', help='Deletion rate (optional)', type = float, default= 1)
    parser_t.add_argument('-m', '--mismatch_rate', help='Mismatch rate (optional)', type = float, default= 1)
    parser_t.add_argument('-max', '--max_len', help='The maximum length for simulated reads', type=int, default= float("inf"))
    parser_t.add_argument('-min', '--min_len', help='The minimum length for simulated reads', type=int, default= 50)
    parser_t.add_argument('-k', '--KmerBias', help='Determine whether to considert Kmer Bias or not', type = int, default= 0)
    parser_t.add_argument('-s', '--strandness', help='Determine the strandness of the simulated reads. Overrides the value profiled in characterization phase. Should be between 0 and 1.', type=float, default=None)
    parser_t.add_argument('--no_model_ir', help='Consider Intron Retention model from characterization step when simulating reads', action='store_true')
    parser_t.add_argument('--perfect', help='Ignore profiles and simulate perfect reads', action='store_true')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 2:
        if args.mode == "genome":
            parser_g.print_help(sys.stderr)
        elif args.mode == "transcriptome":
            parser_t.print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)

    if args.mode == "genome":
        ref_g = args.ref_g
        model_prefix = args.model_prefix
        out = args.output
        number = args.number
        ins_rate = args.insertion_rate
        del_rate = args.deletion_rate
        mis_rate = args.mismatch_rate
        max_readlength = args.max_len
        min_readlength = args.min_len
        median_readlength = args.median_len
        sd_readlength = args.sd_len
        if args.seed:
            random.seed(int(args.seed))
            np.random.seed(int(args.seed))
        if args.perfect:
            perfect = True
        kmer_bias = args.KmerBias
        strandness = args.strandness
        dna_type = args.dna_type

        if dna_type not in ['circular', 'linear']:
            print("\nPlease input proper dna type: linear OR circular\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if kmer_bias < 0:
            print("\nPlease input proper kmer_bias value: k >= 0\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if strandness != None and 0 <= strandness <= 1:
            print("\nPlease input proper strandness value between 0 and 1\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if (median_readlength and not sd_readlength) or (sd_readlength and not median_readlength):
            sys.stderr.write("\nPlease provide both mean and standard deviation of read length!\n\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if max_readlength < min_readlength:
            sys.stderr.write("\nMaximum read length must be longer than Minimum read length!\n\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("ref_g", ref_g)
        print("model_prefix", model_prefix)
        print("out", out)
        print("number", number)
        print("perfect", perfect)
        print("kmer_bias", kmer_bias)
        print("dna_type", dna_type)
        print("strandness", strandness)
        print("sd_readlength", sd_readlength)
        print("median_readlength", median_readlength)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        basename = os.path.basename(out)
        call("mkdir -p " + dir_name, shell=True)

        read_profile(ref_g, None, number, model_prefix, perfect, args.mode, strandness)

        if median_readlength and sd_readlength:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Simulating read length with log-normal distribution\n")
            sys.stdout.flush()
            simulation(args.mode, out, dna_type, perfect, kmer_bias, max_readlength, min_readlength, median_readlength,
                       sd_readlength)
        else:
            simulation(args.mode, out, dna_type, perfect, kmer_bias, max_readlength, min_readlength)

    elif args.mode == "transcriptome":
        ref_g = args.ref_g
        ref_t = args.ref_t
        exp = args.exp
        model_prefix = args.model_prefix
        out = args.output
        number = args.number
        ins_rate = args.insertion_rate
        del_rate = args.deletion_rate
        mis_rate = args.mismatch_rate
        max_readlength = args.max_len
        min_readlength = args.min_len
        kmer_bias = args.KmerBias
        strandness = args.strandness
        if args.perfect:
            perfect = True
        if args.no_model_ir:
            model_ir = False
        dna_type = "transcriptome"

        print("\nrunning the code with following parameters:\n")
        print("ref_g", ref_g)
        print("ref_t", ref_t)
        print("exp", exp)
        print("model_prefix", model_prefix)
        print("out", out)
        print("number", number)
        print("perfect", perfect)
        print("kmer_bias", kmer_bias)
        print("model_ir", model_ir)
        print("dna_type", dna_type)
        print("strandness", strandness)

        dir_name = os.path.dirname(out)
        basename = os.path.basename(out)
        call("mkdir -p " + dir_name, shell=True)

        # Generate log file
        # sys.stdout = open(out + ".log", 'w')
        # Record the command typed to log file
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        read_profile(ref_g, ref_t, number, model_prefix, perfect, args.mode, strandness, exp, model_ir)

        simulation(args.mode, out, dna_type, perfect, kmer_bias, max_readlength, min_readlength, None, None, model_ir)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
    sys.stdout.close()

if __name__ == "__main__":
    main()
