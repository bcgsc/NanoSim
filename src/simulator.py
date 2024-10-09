#!/usr/bin/env python
"""
@author: Chen Yang, Saber Hafezqorani, Ka Ming Nip, and Theodora Lo
This script generates simulated Oxford Nanopore reads (genomic, transcriptomic, and metagenomic).
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
import random
import re
import copy
import argparse
import joblib
from time import strftime
from urllib.request import Request, urlopen
from gzip import GzipFile
import numpy as np
import scipy.stats

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
import model_homopolymer_lengths as model_hp_len
import model_base_qualities as model_base_quals
import math

PYTHON_VERSION = sys.version_info
VERSION = "3.2.1"
PROGRAM = "NanoSim"
AUTHOR = "Chen Yang, Saber Hafezqorani, Ka Ming Nip, and Theodora Lo (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca; shafezqorani@bcgsc.ca; kmnip@bcgsc.ca"

BASES = ['A', 'T', 'C', 'G']


def check_print_progress(sequence_index):
    if sequence_index % 10000 == 0:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Number of reads simulated >> " +
                         str(sequence_index + 1) + "\r")
        sys.stdout.flush()


def list_to_range(input_list, min_l):
    l = [min_l]
    l.extend(input_list)
    output_list = []
    for i in xrange(0, len(l) - 1):
        r = (l[i], l[i + 1])
        output_list.append(r)
    return output_list


def make_cdf(dict_exp, dict_len):
    sum_exp = 0
    match_count = 0
    list_value = []
    for item in dict_exp:
        if item in dict_len:
            match_count += 1
            sum_exp += dict_exp[item]
    if match_count == 0:
        sys.stderr.write("Please make sure transcript IDs in the expression profile match with those in reference transcriptome (example: both Ensembl IDs)\n")
        sys.exit(1)
    for item in dict_exp:
        if item in dict_len:
            value = dict_exp[item] / float(sum_exp)
            list_value.append((item, value))

    sorted_value_list = sorted(list_value, key=lambda x: x[1])
    sorted_only_values = [x[1] for x in sorted_value_list]
    list_cdf = np.cumsum(sorted_only_values)
    ranged_cdf_list = list_to_range(list_cdf, 0)

    ecdf_weight_list = list()
    ecdf_length_list = list()
    for cdf_range, txpt in zip(ranged_cdf_list, sorted_value_list):
        ecdf_weight_list.append(abs(cdf_range[1] - cdf_range[0]))
        tname = txpt[0]
        ecdf_length_list.append((tname, dict_len[tname]))

    return ecdf_length_list, ecdf_weight_list


def ref_len_from_structure(input):
    l = 0
    for item in input:
        if item[0] == "exon":
            l += item[-2]
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
    for i in range(0, count):
        p = random.random()
        for key in IR_markov_model[prev_state]:
            if key[0] <= p < key[1]:
                flag = IR_markov_model[prev_state][key]
                if flag == "IR":
                    flag_ir = True
                list_states.append(flag)
                prev_state = flag
                break

    if flag_ir:
        ref_trx_structure_temp = copy.deepcopy(ref_trx_structure)
        j = -1
        for i in xrange(0, len(ref_trx_structure_temp)):
            if ref_trx_structure_temp[i][0] == "intron":
                j += 1
                if list_states[j] == "IR":
                    ref_trx_structure_temp[i] = ("retained_intron",) + ref_trx_structure_temp[i][1:]
    else:
        ref_trx_structure_temp = ref_trx_structure

    return flag_ir, ref_trx_structure_temp


def extract_read_pos(length, ref_len, ref_trx_structure, polya, buffer=10):
    # The aim is to create a genomic interval object
    # example: iv = HTSeq.GenomicInterval( "chr3", 123203, 127245, "+" )
    # buffer: if the extracted read is within 10 base to the reference 3' end, it's considered as reaching to the end

    # find the length before the first retained intron
    len_before = 0
    for item in ref_trx_structure:
        if item[0] == "exon":
            len_before += item[4]
        elif item[0] == "retained_intron":
            break

    # TODO change the random into something truer
    start_pos = random.randint(0, min(ref_len - length, len_before))  # make sure the retained_intron is included

    list_intervals = []
    ir_list = []
    for item in ref_trx_structure:
        if length == 0:
            break
        chrom = item[1]
        if item[0] in ["exon", "retained_intron"]:
            if start_pos < item[4]:
                start = start_pos + item[2]
                if start + length <= item[3]:
                    end = start + length
                else:
                    end = item[3]
                length -= end - start
                start_pos = 0
                iv = HTSeq.GenomicInterval(chrom, start, end, item[5])
                list_intervals.append(iv)
                if item[0] == "retained_intron":
                    ir_list.append((start, end))
            else:
                start_pos -= item[4]

    if polya and end + buffer >= ref_trx_structure[-1][3]:
        retain_polya = True
    else:
        retain_polya = False

    return list_intervals, retain_polya, ir_list


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
    if flatten:
        return tmp_list.flatten()
    
    return tmp_list


def read_profile(ref_g, number_list, model_prefix, per, mode, strandness, ref_t=None, dna_type=None, abun=None,
                 polya=None, exp=None, model_ir=False, chimeric=False, homopolymer=False, fastq=False):
    # Note var number_list (list) used to be number (int)
    global number_aligned_l, number_unaligned_l, number_segment_list
    global match_ht_list, error_par, trans_error_pr, match_markov_model
    global kde_aligned, kde_ht, kde_ht_ratio, kde_unaligned, kde_aligned_2d
    global seq_dict, seq_len, max_chrom
    global strandness_rate

    if mode == "genome":
        global genome_len
        ref = ref_g
    elif mode == "metagenome":
        global multi_dict_abun, dict_dna_type
        ref = {}
        with open(ref_g, 'r') as genome_list:
            list = genome_list.readlines()
            for genome in list:
                fields = genome.split("\t")
                species = fields[0]
                species = '_'.join(species.split())
                genome = fields[1].strip("\n")
                ref[species] = genome
    else:
        global dict_exp, ecdf_length_list, ecdf_weight_list
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
    dict_dna_type = {}

    # Read in the reference genome/transcriptome/metagenome
    if mode == "metagenome":
        max_chrom = {}
        for species in ref.keys():
            fq_path = ref[species]
            seq_dict[species] = {}
            seq_len[species] = {}
            dict_dna_type[species] = {}
            max_chrom[species] = 0
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in " + species + "\n")
            sys.stdout.flush()

            if fq_path.startswith(("ftp", "http")):
                http_addr = fq_path.replace("ftp://", "http://")
                dir_name = fq_path.strip().split('/')[-1]
                http_complete = http_addr.strip() + '/' + dir_name + '_genomic.fna.gz'
                url_req = Request(http_complete)
                url_req.add_header('Accept-Encoding', 'gzip')
                response = urlopen(url_req)
                with GzipFile(fileobj=response) as f:
                    for line in f:
                        line = str(line, 'utf-8').strip()
                        if line[0] == '>':
                            info = re.split(r'[_\s]\s*', line)
                            chr_name = "-".join(info[1:])
                            seq_dict[species][chr_name.split(".")[0]] = ''
                            seq_len[species][chr_name.split(".")[0]] = 0
                            dict_dna_type[species][chr_name.split(".")[0]] = "linear"  # linear as default
                        else:
                            seq_dict[species][chr_name.split(".")[0]] += line
                            seq_len[species][chr_name.split(".")[0]] += len(line)
                            if seq_len[species][chr_name.split(".")[0]] > max_chrom[species]:
                                max_chrom[species] = seq_len[species][chr_name.split(".")[0]]
            else:
                with open(fq_path, 'r') as infile:
                    for seqN, seqS, seqQ in readfq(infile):
                        info = re.split(r'[_\s]\s*', seqN)
                        chr_name = "-".join(info)
                        seq_dict[species][chr_name.split(".")[0]] = seqS
                        seq_len[species][chr_name.split(".")[0]] = len(seqS)
                        dict_dna_type[species][chr_name.split(".")[0]] = "circular"  # circular as default
                        if len(seqS) > max_chrom[species]:
                            max_chrom[species] = len(seqS)

        if dna_type:  # dna_type is not required when streaming reference genome from RefSeq
            with open(dna_type, 'r') as dna_type_list:
                for line in dna_type_list.readlines():
                    fields = line.split("\t")
                    species = '_'.join(fields[0].split())
                    chr = re.split(r'[_\s]\s*', fields[1].partition(" ")[0])
                    chr_name = "-".join(chr)
                    type = fields[2].strip("\n")

                    if species not in ref:
                        sys.stderr.write("You didn't provide a reference genome for " + species + '\n')
                        sys.exit(1)
                    dict_dna_type[species][chr_name.split(".")[0]] = type
    else:
        max_chrom = 0
        with open(ref, 'r') as infile:
            for seqN, seqS, seqQ in readfq(infile):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                seq_dict[chr_name.split(".")[0]] = seqS
                seq_len[chr_name.split(".")[0]] = len(seqS)
                if len(seqS) > max_chrom:
                    max_chrom = len(seqS)

    # Special files for each mode
    if mode == "genome":
        genome_len = sum(seq_len.values())
        if len(seq_dict) > 1 and dna_type == "circular":
            sys.stderr.write("Do not choose circular if there is more than one chromosome in the genome!\n")
            sys.exit(1)
    elif mode == "metagenome":
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in abundance profile\n")
        sys.stdout.flush()
        with open(abun, 'r') as abun_file:
            header = abun_file.readline()
            number_list = [int(x) for x in header.strip().split('\t')[1:]]
            sample_size = len(number_list)
            samples = ["sample" + str(x) for x in range(sample_size)]
            multi_dict_abun = {sample: {} for sample in samples}
            for line in abun_file.readlines():
                fields = line.split("\t")
                if sample_size != len(fields) - 1:  # abundance file is incorrectly formatted
                    sys.stderr.write("Abundance file is incorrectly formatted. Check that each row has the same number "
                                     "of columns\n")
                    sys.exit(1)

                species = '_'.join(fields[0].split())
                if species not in ref:
                    sys.stderr.write("You didn't provide a reference genome for " + species + '\n')
                    sys.exit(1)

                expected = [float(x) for x in fields[1:]]
                for s_idx in xrange(sample_size):
                    multi_dict_abun[samples[s_idx]][species] = expected[s_idx]

    else:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in expression profile\n")
        sys.stdout.flush()
        dict_exp = {}
        with open(exp, 'r') as exp_file:
            header = exp_file.readline()
            for line in exp_file:
                parts = line.split("\t")
                if len(parts) < 3:
                    sys.stderr.write("Expression profile must contain 3 columns: ID, count, TPM \n")
                    sys.exit(1)
                transcript_id = parts[0].split(".")[0]
                tpm = float(parts[2])
                if tpm > 0:
                    dict_exp[transcript_id] = tpm
        if len(dict_exp) == 0:
            sys.stderr.write("Expression profile contains no TPM values > 0\n")
            sys.exit(1)
        # create the ecdf dict considering the expression profiles
        ecdf_length_list, ecdf_weight_list = make_cdf(dict_exp, seq_len)

        if model_ir:
            global genome_fai, IR_markov_model, dict_ref_structure
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in reference genome and create .fai index file\n")
            sys.stdout.flush()
            # create and read the .fai file of the reference genome
            genome_fai = pysam.Fastafile(ref_g)

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in IR markov model\n")
            sys.stdout.flush()

            IR_markov_model = {}
            with open(model_prefix + "_IR_markov_model", "r") as IR_markov:
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

                    # remove "chr" from chromosome names to be consistent
                    if "chr" in feature.iv.chrom:
                        feature.iv.chrom = feature.iv.chrom.strip("chr")

                    dict_ref_structure[feature_id].append((feature.type, feature.iv.chrom, feature.iv.start,
                                                           feature.iv.end, feature.iv.length, feature.iv.strand))

        if polya:
            global trx_with_polya
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read in list of transcripts with polyA tails\n")
            sys.stdout.flush()
            trx_with_polya = {}
            with open(polya, "r") as trx_list:
                for line in trx_list.readlines():
                    transcript_id = line.strip().split(".")[0]
                    trx_with_polya[transcript_id] = 0

    if per:  # if parameter perfect is used, all reads should be aligned, number_aligned equals total number of reads
        number_aligned_l = number_list
        number_unaligned_l = [0] * len(number_list)
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

        # Read homopolymer length model parameters
        if homopolymer:
            global pw_hp_len, lr_hp_len, hp_mis_rate
            with open(model_prefix + "_hp_lengths_model_parameters.tsv") as homopolymer_length_params:
                pw_hp_len = {}
                lr_hp_len = {}

                # Deal with first two lines of file: mismatch rate and header
                hp_mis_rate = float(re.search("\d+\.?\d*", next(homopolymer_length_params))[0])
                header = next(homopolymer_length_params)
                col_names = header.strip().split("\t")

                # Get parameter values from rest of file
                for line in homopolymer_length_params:
                    fields = line.strip().split("\t")
                    base = fields[0]
                    pw_hp_len[base] = {}
                    lr_hp_len[base] = {}
                    # Iterate through col_names in case number of breakpoints changes and thus, number of
                    # piecewise parameters, changes in the future
                    for i, col_name in enumerate(col_names):
                        if i == 0:
                            continue
                        elif col_name in ["intercept", "slope"]:
                            lr_hp_len[base][col_name] = float(fields[i])
                        else:
                            pw_hp_len[base][col_name] = float(fields[i])

        # Read length of unaligned reads
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read KDF of unaligned reads\n")
        sys.stdout.flush()

        with open(model_prefix + "_reads_alignment_rate", 'r') as u_profile:
            new = u_profile.readline().strip()
            rate = new.split('\t')[1]
            if rate == "100%":
                number_aligned_l = number_list
            else:
                number_aligned_l = [int(round(x * float(rate) / (float(rate) + 1))) for x in number_list]
            number_unaligned_l = [x - y for x, y in zip(number_list, number_aligned_l)]

        if min(number_unaligned_l) > 0:
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
        if mode == "transcriptome":
            kde_aligned_2d = joblib.load(model_prefix + "_aligned_region_2d.pkl")
    else:
        if mode == "transcriptome":
            kde_aligned_2d = joblib.load(model_prefix + "_aligned_region_2d.pkl")
        else:  # genome/metagenome
            kde_aligned = joblib.load(model_prefix + "_aligned_region.pkl")

    # Read chimeric reads information
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read chimeric simulation information\n")
    if chimeric:
        global abun_inflation, kde_gap, segment_mean
        with open(model_prefix + "_chimeric_info") as chimeric_info:
            segment_mean = float(chimeric_info.readline().split('\t')[1])
            if mode == "metagenome":
                abun_inflation = float(chimeric_info.readline().split('\t')[1])
        kde_gap = joblib.load(model_prefix + "_gap_length.pkl")

    # Read base quality model parameters
    if fastq:
        global lognorm_base_qual
        with open(model_prefix + "_base_qualities_model_parameters.tsv") as base_quality_params:
            next(base_quality_params)  # skip header
            lognorm_base_qual = {}
            for line in base_quality_params:
                fields = line.split("\t")
                type = fields[0]
                sd = float(fields[1])
                loc = float(fields[2])
                mu = float(fields[3])
                lognorm_base_qual[type] = {"sd": sd, "loc": loc, "mu": mu}


def add_abundance_var(expected_abun, total_len, var_low, var_high):
    # Order species according to genome size and assign highest variation to species
    # w/ largest genome and lowest variation to species w/ smallest genome
    abun_var = []
    for i in range(len(total_len)):  # Generate % var for each species
        abun_var.append(random.uniform(var_low, var_high))

    abun_var_per_species_dict = {}
    for var, species in zip(sorted(abun_var, key=abs), sorted(total_len, key=lambda k: total_len[k])):
        abun_var_per_species_dict[species] = var

    # Add variation to expected abundances
    abun_with_var = {}
    for species, expected in expected_abun.items():
        abun_with_var[species] = expected + expected * abun_var_per_species_dict[species]

    # Renormalize abundances
    total = sum(abun_with_var.values())
    for species, abun in abun_with_var.items():
        abun_with_var[species] = abun * 100 / total

    return abun_with_var


def mutate_homo(seq, base_quals, k):
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

        hp_length_hist[length][base] += 1

    # Obtaining samples from normal distributions
    for length in hp_length_hist.keys():
        hp_samples[length] = {}
        a_mu, a_sigma, t_mu, t_sigma, c_mu, c_sigma, g_mu, g_sigma = model_hp_len.get_nd_par(length, pw_hp_len, lr_hp_len)

        if hp_length_hist[length]["A"] > 0:
            hp_samples[length]["A"] = np.random.normal(a_mu, a_sigma, hp_length_hist[length]["A"])
        if hp_length_hist[length]["T"] > 0:
            hp_samples[length]["T"] = np.random.normal(t_mu, t_sigma, hp_length_hist[length]["T"])
        if hp_length_hist[length]["C"] > 0:
            hp_samples[length]["C"] = np.random.normal(c_mu, c_sigma, hp_length_hist[length]["C"])
        if hp_length_hist[length]["G"] > 0:
            hp_samples[length]["G"] = np.random.normal(g_mu, g_sigma, hp_length_hist[length]["G"])

    for length in hp_samples.keys():
        for base in hp_samples[length].keys():
            hp_samples[length][base] = [0 if x < 0 else x for x in hp_samples[length][base]]

    # Mutating homopolymers in given sequence
    last_pos = 0
    mutated_seq = ""
    total_hp_size_change = 0
    for hp_info in hp_arr:
        base = hp_info[0]
        ref_hp_start = hp_info[1]
        ref_hp_end = hp_info[2]

        size = int(round(hp_samples[ref_hp_end - ref_hp_start][base][-1]))
        hp_samples[ref_hp_end - ref_hp_start][base] = hp_samples[ref_hp_end - ref_hp_start][base][:-1]

        mutated_hp_with_mis = ""
        mis_pos = []

        for i in xrange(size):
            p = random.random()
            if 0 < p <= hp_mis_rate:
                tmp_bases = list(BASES)
                while True:
                    new_base = random.choice(tmp_bases)
                    if new_base != base:
                        break
                mutated_hp_with_mis += new_base
                mis_pos.append(i)
            else:
                mutated_hp_with_mis += base
            i += 1
        mutated_seq = mutated_seq + seq[last_pos: ref_hp_start] + mutated_hp_with_mis

        if len(base_quals) != 0:  # fastq
            diff = size - (ref_hp_end - ref_hp_start)
            if diff < 0:  # del, remove quals
                for i in xrange(abs(diff)):
                    base_quals.pop(ref_hp_start + total_hp_size_change)

            elif diff > 0:  # ins, add quals
                ins_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["ins"]["sd"], lognorm_base_qual["ins"]["loc"], np.exp(lognorm_base_qual["ins"]["mu"]), diff)
                base_quals = base_quals[:ref_hp_end + total_hp_size_change] + ins_quals + \
                             base_quals[ref_hp_end + total_hp_size_change:]

            if len(mis_pos) != 0:  # mis, change match quals to mis quals
                mis_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["mis"]["sd"], lognorm_base_qual["mis"]["loc"], np.exp(lognorm_base_qual["mis"]["mu"]), 1)
                for i in zip(mis_pos, mis_quals):
                    base_quals[ref_hp_start + total_hp_size_change + i[0]] = i[1]

        total_hp_size_change += size - (ref_hp_end - ref_hp_start)
        last_pos = ref_hp_end

    return mutated_seq + seq[last_pos:], base_quals


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


def assign_species(length_list, seg_list, current_species_base_dict):
    # Deal with chimeric reads first
    seg_list_sorted = sorted(seg_list, reverse=True)
    segs_chimera = sum([x for x in seg_list if x > 1])

    # Sort lengths for non-chimeras to fit in species quota better
    length_list_nonchimera = length_list[segs_chimera:]
    length_list_sorted = length_list[:segs_chimera] + sorted(length_list_nonchimera, reverse=True)

    species_list = [''] * len(length_list)
    bases_to_add = sum(length_list)
    current_bases = sum(current_species_base_dict.values())
    total_bases = bases_to_add + current_bases

    base_quota = {}
    total_abun = sum(dict_abun.values())
    for species, abun in dict_abun.items():
        base_quota[species] = total_bases * abun / total_abun - current_species_base_dict[species]

    length_list_pointer = 0
    pre_species = ''
    num_reads = len(length_list_sorted)
    for seg in seg_list_sorted:
        if length_list_pointer + seg > num_reads:
            break
        for each_seg in range(seg):
            if each_seg == 0:
                available_species = [s for s, q in base_quota.items()
                                     if q - length_list_sorted[length_list_pointer] > 0]
                if len(available_species) == 0:
                    available_species = [s for s, q in base_quota.items() if q > 0]
                species = random.choice(available_species)
            else:
                available_species = [s for s, q in base_quota.items()
                                     if q - length_list_sorted[length_list_pointer] > 0 and s != pre_species]
                p = random.uniform(0, 100)
                if p <= dict_abun_inflated[pre_species] and base_quota[pre_species] > 0:
                    species = pre_species
                elif p > dict_abun_inflated[pre_species] and len(available_species) > 0:
                    species = random.choice(available_species)
                else:
                    available_species = [s for s, q in base_quota.items()
                                         if q - length_list_sorted[length_list_pointer] > 0]
                    if len(available_species) == 0:
                        available_species = [s for s, q in base_quota.items() if q > 0]
                    species = random.choice(available_species)

            species_list[length_list_pointer] = species
            base_quota[species] -= length_list_sorted[length_list_pointer]
            length_list_pointer += 1
            pre_species = species

    return species_list[:length_list_pointer], length_list_sorted[:length_list_pointer], \
           np.array(seg_list_sorted[:length_list_pointer])


def simulation_aligned_metagenome(min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias, fastq, num_simulate,
                                  per=False, chimeric=False):
    # Simulate aligned reads
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    id_begin = '@' if fastq else '>'

    remaining_reads = num_simulate
    if chimeric:
        num_segment = np.random.geometric(1 / segment_mean, num_simulate)
    else:
        num_segment = np.ones(num_simulate, dtype=int)
    remaining_segments = num_segment
    remaining_gaps = remaining_segments - 1
    passed = 0
    current_species_bases = {species: 0 for species in dict_abun.keys()}
    while remaining_reads > 0:
        if per:
            ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments)) if median_l is None else \
                np.random.lognormal(np.log(median_l), sd_l, remaining_segments)
            ref_lengths = [x for x in ref_lengths if min_l <= x <= max_l]
            if len(ref_lengths) == 0:
                continue
        else:
            remainder_lengths = get_length_kde(kde_ht, int(remaining_reads * 1.3), True)
            remainder_lengths = [x for x in remainder_lengths if x >= 0]
            head_vs_ht_ratio_list = get_length_kde(kde_ht_ratio, int(remaining_reads * 1.5))
            head_vs_ht_ratio_list = [x for x in head_vs_ht_ratio_list if 0 <= x <= 1]
            if median_l is None:
                ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments))
            else:
                total_lengths = np.random.lognormal(np.log(median_l + sd_l ** 2 / 2), sd_l, remaining_reads)
                num_current_loop = min(remaining_reads, len(remainder_lengths), len(head_vs_ht_ratio_list))
                ref_lengths = total_lengths[:num_current_loop] - remainder_lengths[:num_current_loop]
            ref_lengths = [x for x in ref_lengths if 0 < x <= max_l]
            if len(ref_lengths) == 0:
                continue

        gap_lengths = get_length_kde(kde_gap, sum(remaining_gaps), True) if sum(remaining_gaps) > 0 else []
        gap_lengths = [max(0, int(x)) for x in gap_lengths]

        # Select strain/species to simulate
        species_pool, ref_lengths, remaining_segments = \
            assign_species(ref_lengths, remaining_segments, current_species_bases)

        is_reversed = random.random() > strandness_rate

        seg_pointer = 0
        gap_pointer = 0
        species_pointer = 0
        for each_read in xrange(len(remaining_segments)):
            segments = remaining_segments[each_read]
            # In case too many ref length was filtered previously
            if (not per and each_read >= min(len(head_vs_ht_ratio_list), len(remainder_lengths))) or \
                seg_pointer + segments > len(ref_lengths):
                break
            ref_length_list = [int(round(ref_lengths[seg_pointer + x])) for x in range(segments)]
            gap_length_list = [int(round(gap_lengths[gap_pointer + x])) for x in range(segments - 1)]
            species_list = [species_pool[species_pointer + x] for x in range(segments)]

            if per:
                seg_pointer += 1
                gap_pointer += 1
                species_pointer += 1
                with total_simulated.get_lock():
                    sequence_index = total_simulated.value
                    total_simulated.value += 1

                # Extract middle region from reference genome
                new_read = ""
                new_read_name = ""
                base_quals = []
                for seg_idx in range(len(ref_length_list)):
                    new_seg, new_seg_name = extract_read("metagenome", ref_length_list[seg_idx], species_list[seg_idx])
                    new_read += new_seg
                    new_read_name += new_seg_name
                    if fastq:
                        base_quals.extend(model_base_quals.predict_base_qualities(lognorm_base_qual["match"]["sd"], lognorm_base_qual["match"]["loc"], np.exp(lognorm_base_qual["match"]["mu"]), ref_length_list[seg_idx]))

                new_read_name = new_read_name + "_perfect_" + str(sequence_index)
                read_mutated = case_convert(new_read)  # not mutated actually, just to be consistent with per == False
                
                if len(read_mutated) < min_l or len(read_mutated) > max_l:
                    continue

                head = 0
                tail = 0

                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_0_" + str(sum(ref_length_list)) + "_0"

            else:
                gap_list = []
                gap_base_qual_list = []
                seg_length_list = []
                seg_error_dict_list = []
                seg_error_count_list = []
                remainder = int(round(remainder_lengths[each_read]))
                head_vs_ht_ratio = head_vs_ht_ratio_list[each_read]

                total = remainder
                restart = False
                for each_ref in ref_length_list:
                    middle, middle_ref, error_dict, error_count = \
                        error_list(each_ref, match_markov_model, match_ht_list, error_par, trans_error_pr, fastq)
                    if total + middle_ref > max_l:
                       restart = True
                       break
                    total += middle_ref
                    seg_length_list.append(middle_ref)
                    seg_error_dict_list.append(error_dict)
                    seg_error_count_list.append(error_count)

                if restart:
                    continue

                for each_gap in gap_length_list:
                    mutated_gap, gap_base_quals = simulation_gap(each_gap, "metagenome", fastq)
                    gap_list.append(mutated_gap)
                    gap_base_qual_list.append(gap_base_quals)
                    gap_length = len(mutated_gap)
                    if total + gap_length > max_l:
                       restart = True
                       break
                    total += gap_length

                if restart or total < min_l or total > max_l:
                    continue

                seg_pointer += segments
                gap_pointer += segments - 1
                species_pointer += segments

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
                num_seg = len(seg_length_list)
                new_seg_list = [None] * num_seg
                read_name_components = list()
                for seg_idx in range(num_seg):
                    new_seg_list[seg_idx], new_seg_name = extract_read("metagenome", seg_length_list[seg_idx], species_list[seg_idx])
                    read_name_components.append(new_seg_name)
                    if seg_idx < len(gap_list):
                        read_name_components.append("gap_" + str(len(gap_list[seg_idx])))
                new_read_name = ';'.join(read_name_components) + "_aligned_" + str(sequence_index)
                
                if num_seg > 1:
                    new_read_name += "_chimeric"
                
                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_" + str(head) + \
                                 "_" + ";".join(str(x) for x in seg_length_list) + \
                                 "_" + str(tail)
                    
                read_mutated = ""
                base_quals = []
                for seg_idx in range(num_seg):
                    # Mutate read
                    new_seg = case_convert(new_seg_list[seg_idx])
                    seg_mutated, seg_base_quals = \
                        mutate_read(new_seg, new_read_name, out_error, seg_error_dict_list[seg_idx],
                                    seg_error_count_list[seg_idx], fastq, kmer_bias)

                    if kmer_bias:
                        seg_mutated, seg_base_quals = mutate_homo(seg_mutated, seg_base_quals, kmer_bias)
                    read_mutated += seg_mutated
                    base_quals.extend(seg_base_quals)
                    if seg_idx < len(gap_list):
                        read_mutated += gap_list[seg_idx]
                        base_quals.extend(gap_base_qual_list[seg_idx])

                    # Update base level abundance info
                    current_species_bases[species_list[seg_idx]] += len(new_seg)

                if fastq:  # Get head/tail qualities and add to base_quals
                    ht_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["ht"]["sd"], lognorm_base_qual["ht"]["loc"], np.exp(lognorm_base_qual["ht"]["mu"]), head + tail)
                    base_quals = ht_quals[:head] + base_quals + ht_quals[head:]

            # Add head and tail region
            read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                           ''.join(np.random.choice(BASES, tail))

            if len(read_mutated) < min_l or len(read_mutated) > max_l:
                continue

            # Reverse complement half of the reads
            if is_reversed:
                read_mutated = reverse_complement(read_mutated)
                base_quals.reverse()

            out_reads.write(id_begin + new_read_name + '\n')
            out_reads.write(read_mutated + '\n')

            if fastq:
                out_reads.write("+\n")
                out_quals = "".join([chr(qual + 33) for qual in base_quals])
                out_reads.write(out_quals + "\n")

            check_print_progress(sequence_index)

            passed += 1

        remaining_reads = num_simulate - passed
        remaining_segments = num_segment[passed:]
        remaining_gaps = remaining_segments - 1

    sys.stdout.write('\n')
    out_reads.close()
    out_error.close()


def simulation_aligned_transcriptome(model_ir, out_reads, out_error, kmer_bias, basecaller, num_simulate, polya, fastq,
                                     per=False, uracil=False):

    if basecaller == "albacore":
        polya_len_dist_scale = 2.409858743694814
    else:
        polya_len_dist_scale = 4.168299657168961
    
    # inner function to return a polyA tail length    
    def get_polya_len():
        return int(scipy.stats.expon.rvs(loc=2.0, scale=polya_len_dist_scale))
                                         
    # Simulate aligned reads
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    if fastq:
        id_begin = "@"
    else:
        id_begin = ">"

    if model_ir:
        # check whether chrom names contains "chr" or not.
        flag_chrom = False
        for item in genome_fai.references:
            if "chr" in item:
                flag_chrom = True
                break

    remainder_l = get_length_kde(kde_ht, num_simulate, True)
    head_vs_ht_ratio_temp = get_length_kde(kde_ht_ratio, num_simulate)
    head_vs_ht_ratio_l = [1 if x > 1 else x for x in head_vs_ht_ratio_temp]
    head_vs_ht_ratio_l = [0 if x < 0 else x for x in head_vs_ht_ratio_l]

    simulated = 0
    sampled_2d_lengths = get_length_kde(kde_aligned_2d, num_simulate, False, False) # initial sample from the KDE
    trx_sampled = set() # track transcript IDs that are associated with the KDE sample
    
    while simulated < num_simulate:
        while True:            
            # select a random reference transcript
            ref_trx, ref_trx_len = random.choices(ecdf_length_list, weights=ecdf_weight_list, k=1)[0]
            
            # if this transcript was previously associated with the currect KDE sample
            if ref_trx in trx_sampled:
                # draw a new sample from the KDE
                sampled_2d_lengths = get_length_kde(kde_aligned_2d, num_simulate, False, False)
                
                # track a new set of transcript IDs
                trx_sampled = set()

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
        
        # associate the transcript ID to the KDE sample
        trx_sampled.add(ref_trx)
               
        trx_has_polya = polya and ref_trx in trx_with_polya     
        is_reversed = random.random() > strandness_rate
            
        if per:
            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            new_read, ref_start_pos, retain_polya = extract_read_trx(ref_trx, ref_len_aligned, trx_has_polya)
            new_read_name = ref_trx + "_" + str(ref_start_pos) + "_perfect_" + str(sequence_index)
            read_mutated = case_convert(new_read)  # not mutated actually, just to be consistent with per == False

            if fastq:
                base_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["match"]["sd"], lognorm_base_qual["match"]["loc"], np.exp(lognorm_base_qual["match"]["mu"]), ref_len_aligned)
            else:
                base_quals = []

            head = 0
            tail = 0
            
            if is_reversed:
                new_read_name += "_R"
            else:
                new_read_name += "_F"
            
            if retain_polya:
                polya_len = get_polya_len()
                if polya_len > 0:
                    read_mutated += "A" * polya_len
            else:
                polya_len = 0
                        
            new_read_name += "_0_" + str(ref_len_aligned) + "_" + str(polya_len)
            
        else:
            middle_read, middle_ref, error_dict, error_count = error_list(ref_len_aligned, match_markov_model,
                                                                          match_ht_list, error_par, trans_error_pr,
                                                                          fastq)

            if middle_ref > ref_trx_len:
                continue

            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            ir_list = []
            if model_ir:
                ir_flag, ref_trx_structure_new = update_structure(dict_ref_structure[ref_trx], IR_markov_model)
                if ir_flag:
                    list_iv, retain_polya, ir_list = extract_read_pos(middle_ref, ref_trx_len, ref_trx_structure_new,
                                                                      trx_has_polya)
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
                        new_read += genome_fai.fetch(chrom, start, end)  # len(new_read) > middle_ref
                    if flag:
                        continue
                    ref_start_pos = list_iv[0].start
 
                    if interval.strand == '-':  # Keep the read direction the same as reference transcripts
                        new_read = reverse_complement(new_read)

                    if fastq:  # since len(new_read) > middle_ref if IR, add more match quals for retained intron
                        error_count["match"] += len(new_read) - middle_ref
                else:
                    new_read, ref_start_pos, retain_polya = extract_read_trx(ref_trx, middle_ref, trx_has_polya)

            else:
                new_read, ref_start_pos, retain_polya = extract_read_trx(ref_trx, middle_ref, trx_has_polya)

            new_read_name = str(ref_trx) + "_" + str(ref_start_pos) + "_aligned_" + str(sequence_index)
            if len(ir_list) > 0:
                new_read_name += "_RetainedIntron_"
                for ir_tuple in ir_list:
                    new_read_name += '-'.join(str(x) for x in ir_tuple) + ';'
            
            if is_reversed:
                new_read_name += "_R"
            else:
                new_read_name += "_F"
            
            # start HD len simulation
            remainder = int(remainder_l[simulated])
            head_vs_ht_ratio = head_vs_ht_ratio_l[simulated]
            
            if remainder == 0:
                head = 0
                tail = 0
            else:
                head = int(round(remainder * head_vs_ht_ratio))
                tail = remainder - head
            # end HD len simulation
            
            if retain_polya:
                polya_len = get_polya_len()
            else:
                polya_len = 0
            
            
            new_read_name += "_" + str(head) + \
                             "_" + str(middle_ref) + \
                             "_" + str(tail + polya_len)
            
            # Mutate read
            new_read = case_convert(new_read)
            read_mutated, base_quals = mutate_read(new_read, new_read_name, out_error, error_dict, error_count,
                                                   fastq, kmer_bias)
            if kmer_bias:
                read_mutated, base_quals = mutate_homo(read_mutated, base_quals, kmer_bias)
            
            if polya_len > 0:
                read_mutated += "A" * polya_len

        if fastq:  # Get head/tail qualities
            ht_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["ht"]["sd"], lognorm_base_qual["ht"]["loc"], np.exp(lognorm_base_qual["ht"]["mu"]), head + tail + polya_len)
            for a in xrange(polya_len):
                base_quals.append(ht_quals.pop())
            base_quals = ht_quals[:head] + base_quals + ht_quals[head:]

        # Add head and tail region
        read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + ''.join(np.random.choice(BASES, tail))
        
        # Reverse complement according to strandness rate
        if is_reversed:
            read_mutated = reverse_complement(read_mutated)
            base_quals.reverse()

        out_reads.write(id_begin + new_read_name + '\n')

        if uracil:
            read_mutated = read_mutated.translate(trantab)

        out_reads.write(read_mutated + '\n')

        if fastq:
            out_reads.write("+\n")
            out_quals = "".join([chr(qual + 33) for qual in base_quals])
            out_reads.write(out_quals + "\n")

        check_print_progress(sequence_index)

        simulated += 1

    sys.stdout.write('\n')
    out_reads.close()
    out_error.close()


def simulation_aligned_genome(dna_type, min_l, max_l, median_l, sd_l, out_reads, out_error, kmer_bias, fastq,
                              num_simulate, per=False, chimeric=False):

    # Simulate aligned reads
    out_reads = open(out_reads, "w")
    out_error = open(out_error, "w")

    id_begin = '@' if fastq else '>'

    remaining_reads = num_simulate
    if chimeric:
        num_segment = np.random.geometric(1/segment_mean, num_simulate)
    else:
        num_segment = np.ones(num_simulate, dtype=int)
    remaining_segments = num_segment
    remaining_gaps = remaining_segments - 1
    passed = 0
    while remaining_reads > 0:
        if per:
            ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments)) if median_l is None else \
                np.random.lognormal(np.log(median_l), sd_l, remaining_segments)
            ref_lengths = [x for x in ref_lengths if min_l <= x <= max_l]
        else:
            remainder_lengths, head_vs_ht_ratio_list = get_lengths_and_ht_ratios(remaining_reads)
            if median_l is None:
                ref_lengths = get_length_kde(kde_aligned, sum(remaining_segments))
            else:
                total_lengths = np.random.lognormal(np.log(median_l + sd_l ** 2 / 2), sd_l, remaining_reads)
                num_current_loop = min(remaining_reads, len(remainder_lengths), len(head_vs_ht_ratio_list))
                ref_lengths = total_lengths[:num_current_loop] - remainder_lengths[:num_current_loop]
            ref_lengths = [x for x in ref_lengths if 0 < x <= max_l]

        gap_lengths = get_length_kde(kde_gap, sum(remaining_gaps), True) if sum(remaining_gaps) > 0 else []
        gap_lengths = [max(0, int(x)) for x in gap_lengths]

        seg_pointer = 0
        gap_pointer = 0
        for each_read in xrange(remaining_reads):
            # check if the total length fits the criteria
            segments = remaining_segments[each_read]
            # In case too many ref length was filtered previously
            if seg_pointer + segments > len(ref_lengths):
                break
            ref_length_list = [int(ref_lengths[seg_pointer + x]) for x in range(segments)]
            gap_length_list = [int(gap_lengths[gap_pointer + x]) for x in range(segments - 1)]

            is_reversed = random.random() > strandness_rate

            if per:
                seg_pointer += 1
                gap_pointer += 1
                with total_simulated.get_lock():
                    sequence_index = total_simulated.value
                    total_simulated.value += 1

                # Extract middle region from reference genome
                new_read = ""
                new_read_name = ""
                base_quals = []
                for each_ref in ref_length_list:
                    new_seg, new_seg_name = extract_read(dna_type, each_ref)
                    new_read += new_seg
                    new_read_name += new_seg_name
                    if fastq:
                        base_quals.extend(model_base_quals.predict_base_qualities(lognorm_base_qual["match"]["sd"], lognorm_base_qual["match"]["loc"], np.exp(lognorm_base_qual["match"]["mu"]), each_ref))

                new_read_name = new_read_name + "_perfect_" + str(sequence_index)
                read_mutated = case_convert(new_read)  # not mutated actually, just to be consistent with per == False

                head = 0
                tail = 0
                
                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_0_" + str(sum(ref_length_list)) + "_0"
                
            else:
                gap_list = []
                gap_base_qual_list = []
                seg_length_list = []
                seg_error_dict_list = []
                seg_error_count_list = []
                remainder = int(remainder_lengths[each_read])
                head_vs_ht_ratio = head_vs_ht_ratio_list[each_read]

                total = remainder
                for each_gap in gap_length_list:
                    mutated_gap, gap_base_quals = simulation_gap(each_gap, dna_type, fastq)
                    gap_list.append(mutated_gap)
                    gap_base_qual_list.append(gap_base_quals)
                for each_ref in ref_length_list:
                    middle, middle_ref, error_dict, error_count = \
                        error_list(each_ref, match_markov_model, match_ht_list, error_par, trans_error_pr, fastq)
                    total += middle
                    seg_length_list.append(middle_ref)
                    seg_error_dict_list.append(error_dict)
                    seg_error_count_list.append(error_count)

                if total < min_l or total > max_l:
                    continue

                seg_pointer += segments
                gap_pointer += segments - 1

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
                num_seg = len(seg_length_list)
                new_seg_list = [None] * num_seg
                read_name_components = [None] * num_seg
                for seg_idx in range(num_seg):
                    new_seg_list[seg_idx], read_name_components[seg_idx] = extract_read(dna_type, seg_length_list[seg_idx])
                new_read_name = ';'.join(read_name_components) + "_aligned_" + str(sequence_index)
                
                if num_seg > 1:
                    new_read_name += "_chimeric"
                
                if is_reversed:
                    new_read_name += "_R"
                else:
                    new_read_name += "_F"
                
                new_read_name += "_" + str(head) + \
                                 "_" + ";".join(str(x) for x in seg_length_list) + \
                                 "_" + str(tail)
                
                read_mutated = ""
                base_quals = []
                for seg_idx in range(num_seg):
                    # Mutate read
                    new_seg = case_convert(new_seg_list[seg_idx])
                    seg_mutated, seg_base_quals = \
                        mutate_read(new_seg, new_read_name, out_error, seg_error_dict_list[seg_idx],
                                    seg_error_count_list[seg_idx], fastq, kmer_bias)

                    if kmer_bias:
                        seg_mutated, seg_base_quals = mutate_homo(seg_mutated, seg_base_quals, kmer_bias)
                    read_mutated += seg_mutated
                    base_quals.extend(seg_base_quals)
                    if seg_idx < len(gap_list):
                        read_mutated += gap_list[seg_idx]
                        base_quals.extend(gap_base_qual_list[seg_idx])

                if fastq:  # Get head/tail qualities and add to base_quals
                    ht_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["ht"]["sd"], lognorm_base_qual["ht"]["loc"], np.exp(lognorm_base_qual["ht"]["mu"]), head + tail)
                    base_quals = ht_quals[:head] + base_quals + ht_quals[head:]

            # Add head and tail region
            read_mutated = ''.join(np.random.choice(BASES, head)) + read_mutated + \
                           ''.join(np.random.choice(BASES, tail))

            if len(read_mutated) < min_l or len(read_mutated) > max_l:
                continue

            # Reverse complement half of the reads
            if is_reversed:
                read_mutated = reverse_complement(read_mutated)
                base_quals.reverse()

            out_reads.write(id_begin + new_read_name + '\n')
            out_reads.write(read_mutated + '\n')

            if fastq:
                out_reads.write("+\n")
                out_quals = "".join([chr(qual + 33) for qual in base_quals])
                out_reads.write(out_quals + "\n")

            check_print_progress(sequence_index)

            passed += 1

        remaining_reads = num_simulate - passed
        remaining_segments = num_segment[passed:]
        remaining_gaps = remaining_segments - 1

    out_reads.close()
    out_error.close()

def get_lengths_and_ht_ratios(remaining_reads):
    "Return the remainer_lengths and head_vs_ht_ratio_list lists, ensuring that the lengths are at least remaining_reads long to avoid IndexErrors"
    remainder_multiplier = 1.3 # Settings previously set by Chen
    ht_ratio_multiplier = 1.5
    
    remainder_lengths, head_vs_ht_ratio_list = [],[]
    
    iteration_count = 0
    
    while len(remainder_lengths) < remaining_reads or len(head_vs_ht_ratio_list) < remaining_reads:
        if iteration_count > 50:
            print("Warning - After 50 iterations, insufficient remainder_lengths and/or head_vs_ht_ratio_list. "
                  "Exiting to avoid infinite loop. Please check the KDE training and/or the min/max length for simulated reads "
                  "(Making these length specifications more permissive may avoid this issue.). ", file=sys.stderr)
            break
        remainder_lengths = get_length_kde(kde_ht, int(remaining_reads * remainder_multiplier), True)
        remainder_lengths = [x for x in remainder_lengths if x >= 0]
        head_vs_ht_ratio_list = get_length_kde(kde_ht_ratio, int(remaining_reads * ht_ratio_multiplier))
        head_vs_ht_ratio_list = [x for x in head_vs_ht_ratio_list if 0 <= x <= 1]
        remainder_multiplier *= 1.5 # Increase the ratios to get longer lists, ensure they are longer than the remaining_reads int
        ht_ratio_multiplier *= 1.5
        iteration_count += 1

    return remainder_lengths,head_vs_ht_ratio_list


def simulation_unaligned(dna_type, min_l, max_l, median_l, sd_l, out_reads, fastq, num_simulate, uracil):
    out_reads = open(out_reads, "w")

    if fastq:
        id_begin = "@"
    else:
        id_begin = ">"

    remaining_reads = num_simulate
    passed = 0
    while remaining_reads > 0:
        # if the median length and sd is set, use log normal distribution for simulation
        ref_l = get_length_kde(kde_unaligned, remaining_reads) if median_l is None else \
            np.random.lognormal(np.log(median_l), sd_l, remaining_reads)

        for j in xrange(len(ref_l)):
            # check if the total length fits the criteria
            ref = int(ref_l[j])

            unaligned, middle_ref, error_dict, error_count = unaligned_error_list(ref, error_par)

            if middle_ref < min_l or middle_ref > max_l:
                continue

            with total_simulated.get_lock():
                sequence_index = total_simulated.value
                total_simulated.value += 1

            new_read, new_read_name = extract_read(dna_type, middle_ref)
            new_read_name = new_read_name + "_unaligned_" + str(sequence_index)
            # Change lowercase to uppercase and replace N with any base
            new_read = case_convert(new_read)
            # no quals returned here since unaligned quals are not based on mis/ins/match qual distributions
            read_mutated, _ = mutate_read(new_read, new_read_name, None, error_dict, error_count, False, False)
            
            if len(read_mutated) < min_l or len(read_mutated) > max_l:
                continue

            if fastq:
                base_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["unmapped"]["sd"], lognorm_base_qual["unmapped"]["loc"], np.exp(lognorm_base_qual["unmapped"]["mu"]), len(read_mutated))
            else:
                base_quals = []

            # Reverse complement some of the reads based on direction information
            p = random.random()
            if p > strandness_rate:
                read_mutated = reverse_complement(read_mutated)
                new_read_name += "_R"
                base_quals.reverse()
            else:
                new_read_name += "_F"

            out_reads.write(id_begin + new_read_name + "_0_" + str(middle_ref) + "_0" + '\n')
            if uracil:
                read_mutated = read_mutated.traslate(trantab)
            out_reads.write(read_mutated + "\n")

            if fastq:
                out_reads.write("+\n")
                out_quals = "".join([chr(qual + 33) for qual in base_quals])
                out_reads.write(out_quals + "\n")

            check_print_progress(sequence_index)
            
            passed += 1

        remaining_reads = num_simulate - passed
    out_reads.close()


def simulation_gap(ref, dna_type, fastq):
    if ref == 0:
        return '', []

    unaligned, middle_ref, error_dict, error_count = unaligned_error_list(ref, error_par)
    new_gap, new_gap_name = extract_read(dna_type, middle_ref)
    new_gap = case_convert(new_gap)

    # no quals returned here since unaligned quals are not based on mis/ins/match qual distributions
    gap_mutated, _ = mutate_read(new_gap, new_gap_name, None, error_dict, error_count, False, False)

    if fastq:
        base_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["unmapped"]["sd"], lognorm_base_qual["unmapped"]["loc"], np.exp(lognorm_base_qual["unmapped"]["mu"]), len(gap_mutated))
    else:
        base_quals = []

    return gap_mutated, base_quals


def simulation(mode, out, dna_type, per, kmer_bias, basecaller, max_l, min_l, num_threads, fastq,
               median_l=None, sd_l=None, model_ir=False, uracil=False, polya=None, chimeric=False):
    global total_simulated  # Keeps track of number of reads that have been simulated so far
    total_simulated = mp.Value("i", 0, lock=True)

    mp.set_start_method('fork', force=True)  # TODO: Remove this later, for testing purposes
    # Start simulation
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of aligned reads\n")
    sys.stdout.flush()
    if fastq:
        ext = ".fastq"
    else:
        ext = ".fasta"

    procs = []
    aligned_subfiles = []
    error_subfiles = []
    num_simulate = int(number_aligned / num_threads)

    for i in range(num_threads):
        np.random.seed()
        random.seed()
        aligned_subfile = out + "_aligned_reads{}".format(i) + ext
        error_subfile = out + "_error_profile{}".format(i)
        aligned_subfiles.append(aligned_subfile)
        error_subfiles.append(error_subfile)
        if i == num_threads - 1:  # Last process will simulate the remaining reads
            num_simulate += number_aligned % num_threads

        if mode == "genome":
            p = mp.Process(target=simulation_aligned_genome,
                           args=(dna_type, min_l, max_l, median_l, sd_l, aligned_subfile, error_subfile,
                                 kmer_bias, fastq, num_simulate, per, chimeric))
            procs.append(p)
            p.start()

        elif mode == "metagenome":
            p = mp.Process(target=simulation_aligned_metagenome,
                           args=(min_l, max_l, median_l, sd_l, aligned_subfile, error_subfile, kmer_bias,
                                 fastq, num_simulate, per, chimeric))
            procs.append(p)
            p.start()

        else:
            p = mp.Process(target=simulation_aligned_transcriptome,
                           args=(model_ir, aligned_subfile, error_subfile, kmer_bias, basecaller, num_simulate, polya,
                                 fastq, per, uracil))
            procs.append(p)
            p.start()

    for p in procs:
        p.join()

    sys.stdout.write('\n')  # Start a new line because the "Number of reads simulated" is not returned
    # Merging aligned reads subfiles and error subfiles
    with open(out + "_aligned_reads" + ext, 'w') as out_aligned_reads:
        for fname in aligned_subfiles:
            with open(fname) as infile:
                out_aligned_reads.write(infile.read())
    for fname in aligned_subfiles:
        os.remove(fname)

    with open(out + "_aligned_error_profile", 'w') as out_error:
        out_error.write("Seq_name\tSeq_pos\terror_type\terror_length\tref_base\tseq_base\n")
        for fname in error_subfiles:
            with open(fname) as infile:
                out_error.write(infile.read())
    for fname in error_subfiles:
        os.remove(fname)

    # Simulate unaligned reads, if per, number_unaligned = 0, taken care of in read_ecdf
    if not per:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Start simulation of random reads\n")
        sys.stdout.flush()
        unaligned_subfiles = []
        # unaligned_error_subfiles = []
        num_simulate = int(number_unaligned / num_threads)
        for i in range(num_threads):
            unaligned_subfile = out + "_unaligned_reads{}".format(i) + ext
            # unaligned_error_subfile = out + "_unaligned_error_profile{}".format(i) + ext
            unaligned_subfiles.append(unaligned_subfile)
            # unaligned_error_subfiles.append(unaligned_error_subfile)
            if i == num_threads - 1:
                num_simulate += number_unaligned % num_threads

            # Dividing number of unaligned reads that need to be simulated amongst the number of processes
            p = mp.Process(target=simulation_unaligned,
                           args=(dna_type, min_l, max_l, median_l, sd_l, unaligned_subfile, fastq, num_simulate, uracil))
            procs.append(p)
            p.start()

        for p in procs:
            p.join()

        sys.stdout.write('\n')  # Start a new line because the "Number of reads simulated" is not returned
        # Merging unaligned reads subfiles and error subfiles
        with open(out + "_unaligned_reads" + ext, 'w') as out_unaligned_reads:
            for fname in unaligned_subfiles:
                with open(fname) as infile:
                    out_unaligned_reads.write(infile.read())
        for fname in unaligned_subfiles:
            os.remove(fname)


def reverse_complement(seq):
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    seq_list = list(seq)
    reverse_seq_list = reversed([comp.get(base, base) for base in seq_list])
    reverse_seq = ''.join(reverse_seq_list)
    return reverse_seq


def extract_read_trx(key, length, trx_has_polya, buffer=10):
    # buffer: if the extracted read is within 10 base to the reference 3' end, it's considered as reaching to the end
    # TODO change the random into something truer
    ref_pos = random.randint(0, seq_len[key] - length)
    new_read = seq_dict[key][ref_pos: ref_pos + length]
    retain_polya = False
    if trx_has_polya and ref_pos + length + buffer >= seq_len[key]:  # Read reaches end of transcript
        retain_polya = True
    return new_read, ref_pos, retain_polya


def extract_read(dna_type, length, s=None):
    if dna_type == "transcriptome":
        while True:
            key = random.choice(list(seq_len.keys()))  # added "list" thing to be compatible with Python v3
            if length < seq_len[key]:
                ref_pos = random.randint(0, seq_len[key] - length)
                new_read = seq_dict[key][ref_pos: ref_pos + length]
                new_read_name = key + "_" + str(ref_pos)
                break
        return new_read, new_read_name
    elif dna_type == "metagenome":
        if not s:
            s = random.choice(list(seq_len.keys()))  # added "list" thing to be compatible with Python v3
        
        key = random.choice(list(seq_len[s].keys()))
        key_seq_len = seq_len[s][key]
        
        if length > key_seq_len:
            # the chromosome selected is too short
            
            # find all chromosomes that are longer
            longer_chroms = list()
            longer_chroms_target = list()
            for tmp_s in seq_len:
                tmp_dict = seq_len[tmp_s]
                for tmp_key in tmp_dict:
                    if length < tmp_dict[tmp_key]:
                        if tmp_s == s:
                            longer_chroms_target.append((tmp_s, tmp_key))
                        else:
                            longer_chroms.append((tmp_s, tmp_key))
            
            assert len(longer_chroms) > 0 or len(longer_chroms_target) > 0 # otherwise there is a problem
            # select a random chromosome from this list
            if longer_chroms_target:
                s, key = random.choice(longer_chroms_target)
            else:
                s, key = random.choice(longer_chroms)
                print("Warning: chosen species/strain is shorter than the simulated read length, randomly selected another species/strain. "
                      "It is recommended to run read_analysis.py quantify after simulation to check abundances in simulated reads.",
                      file = sys.stderr)
            key_seq_len = seq_len[s][key]
        
        if dict_dna_type[s][key] == "circular":
            ref_pos = random.randint(0, key_seq_len)
            if length + ref_pos > key_seq_len:
                new_read = seq_dict[s][key][ref_pos:]
                new_read = new_read + seq_dict[s][key][0: length - key_seq_len + ref_pos]
            else:
                new_read = seq_dict[s][key][ref_pos: ref_pos + length]
        else:
            ref_pos = random.randint(0, key_seq_len - length)
            new_read = seq_dict[s][key][ref_pos: ref_pos + length]
        new_read_name = s + '-' + key + "_" + str(ref_pos)
            
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
            # This is designed for genomes with multiple chromosomes which varies a lot in lengths

            while True:
                new_read = ""
                ref_pos = random.randint(0, genome_len)
                for key in seq_len:
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
    e_count = {"match": 0, "mis": 0, "ins": 0}  # Not used; added to be consistent with error_list()
    if m_ref == 0:
        return l_new, middle_ref, e_dict, e_count
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

    return l_new, middle_ref, e_dict, e_count


def error_list(m_ref, m_model, m_ht_list, error_p, trans_p, fastq):
    # l_old is the original length, and l_new is used to control the new length after introducing errors
    l_new = m_ref
    pos = 0
    e_dict = {}
    middle_ref = m_ref
    prev_error = "start"
    e_count = {"mis": 0, "ins": 0, "match": 0}

    # The first match come from m_ht_list
    p = random.random()
    k1 = list(m_ht_list.keys())[0]
    for k2, v2 in m_ht_list[k1].items():
        if k2[0] < p <= k2[1]:
            prev_match = int(np.floor((p - k2[0]) / (k2[1] - k2[0]) * (v2[1] - v2[0]) + v2[0]))
            if prev_match < 2:
                prev_match = 2
    pos += prev_match
    if fastq:
        if prev_match > middle_ref: 
            e_count["match"] += middle_ref
        else:
            e_count["match"] += prev_match

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

        if fastq:
            if error == "mis" or error == "ins":
                e_count[error] += step

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

        if fastq:
            e_count["match"] += step

        if pos + prev_match > middle_ref:
            l_new += pos + prev_match - middle_ref
            middle_ref = pos + prev_match

        pos += prev_match
        if prev_match == 0:
            prev_error += "0"

    return l_new, middle_ref, e_dict, e_count


def mutate_read(read, read_name, error_log, e_dict, e_count, fastq, k):
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
                if not (hp_end <= err_start or err_end <= hp_start):  # Lands in hp; remove err
                    hp_err = True
                    if fastq:  # Convert err qual to match qual
                        if err != "ins":
                            e_count["match"] += e_dict[err_start][1]
                        if err != "del":
                            e_count[err] -= e_dict[err_start][1]
                    break

            if not hp_err:
                new_e_dict[err_start] = [err, e_dict[err_start][1]]

    else:
        new_e_dict = e_dict

    if fastq:  # Sample base qualities for mis/ins/match
        mis_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["mis"]["sd"], lognorm_base_qual["mis"]["loc"], np.exp(lognorm_base_qual["mis"]["mu"]), e_count["mis"])
        ins_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["ins"]["sd"], lognorm_base_qual["ins"]["loc"], np.exp(lognorm_base_qual["ins"]["mu"]), e_count["ins"])
        match_quals = model_base_quals.predict_base_qualities(lognorm_base_qual["match"]["sd"], lognorm_base_qual["match"]["loc"], np.exp(lognorm_base_qual["match"]["mu"]), e_count["match"])

    # Mutate read
    quals = []
    prev = len(read)
    for key in sorted(new_e_dict.keys(), reverse=True):
        val = new_e_dict[key]
        key = math.ceil(key)  # Ceil instead of round for consistent match calculations during base qual sim
        err_quals = []

        if val[0] == "mis":
            ref_base = read[key: key + val[1]]
            new_bases = ""
            for i in xrange(val[1]):
                tmp_bases = list(BASES)
                tmp_bases.remove(read[key + i])
                # tmp_bases.remove(read[key]) ## Edited this part for testing
                new_base = random.choice(tmp_bases)
                new_bases += new_base
                if fastq:
                    err_quals.append(mis_quals.pop())

            new_read = read[:key] + new_bases + read[key + val[1]:]
            err_end = key + val[1]

        elif val[0] == "del":
            new_bases = val[1] * "-"
            ref_base = read[key: key + val[1]]
            new_read = read[: key] + read[key + val[1]:]
            err_end = key + val[1]

        elif val[0] == "ins":
            ref_base = val[1] * "-"
            new_bases = ""
            for i in xrange(val[1]):
                new_base = random.choice(BASES)
                new_bases += new_base
                if fastq:
                    err_quals.append(ins_quals.pop())
            new_read = read[:key] + new_bases + read[key:]
            err_end = key

        if fastq:
            if err_end != prev:  # Match after error
                for j in xrange(prev - err_end):
                    quals.append(match_quals.pop())
            quals += err_quals

        read = new_read
        prev = key

        if val[0] != "match" and error_log:
            error_log.write(read_name + "\t" + str(key) + "\t" + val[0] + "\t" + str(val[1]) +
                            "\t" + ref_base + "\t" + new_bases + "\n")

    if fastq:  # Add first match quals
        while len(match_quals) > 0:
            quals.append(match_quals.pop())

    quals.reverse()
    return read, quals


def inflate_abun(original_dict, inflated_species):
    rest_abun = (1 - original_dict[inflated_species]) * abun_inflation
    inflated_prob = 1 - rest_abun

    return inflated_prob


def main():
    global number_aligned, number_unaligned
    parser = argparse.ArgumentParser(
        description=dedent('''
        Simulation step
        -----------------------------------------------------------
        Given error profiles, reference genome, metagenome,
        and/or transcriptome, simulate ONT DNA or RNA reads
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='NanoSim ' + VERSION)
    subparsers = parser.add_subparsers(help="You may run the simulator on genome, transcriptome, or metagenome mode.",
                                       dest='mode', description=dedent('''
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
    parser_g.add_argument('-med', '--median_len', help='The median read length (Default = None), Note: this simulation '
                                                       'is not compatible with chimeric reads simulation',
                          type=int, default=None)
    parser_g.add_argument('-sd', '--sd_len', help='The standard deviation of read length in log scale (Default = None),'
                                                  ' Note: this simulation is not compatible with chimeric reads '
                                                  'simulation', type=float, default=None)
    parser_g.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_g.add_argument('-hp', '--homopolymer', help='Simulate homopolymer lengths (Default = False)',
                          action='store_true', default=False)
    parser_g.add_argument('-k', '--KmerBias', help='Minimum homopolymer length to simulate homopolymer contraction and '
                                                   'expansion events in, a typical k is 5',
                          type=int, default=None)
    parser_g.add_argument('-s', '--strandness', help='Proportion of sense sequences. Overrides the value '
                                                      'profiled in characterization stage. Should be between 0 and 1',
                          type=float, default=None)
    parser_g.add_argument('-dna_type', help='Specify the dna type: circular OR linear (Default = linear)',
                          choices=["linear", "circular"], default="linear")
    parser_g.add_argument('--perfect', help='Ignore error profiles and simulate perfect reads', action='store_true',
                          default=False)
    parser_g.add_argument('--fastq', help='Output fastq files instead of fasta files', action='store_true',
                          default=False)
    parser_g.add_argument('--chimeric', help='Simulate chimeric reads', action='store_true', default=False)
    parser_g.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                          default=1)

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode")
    parser_t.add_argument('-rt', '--ref_t', help='Input reference transcriptome', required=True)
    parser_t.add_argument('-rg', '--ref_g', help='Input reference genome, required if intron retention simulation is '
                                                 'on', default='')
    parser_t.add_argument('-e', '--exp', help='Expression profile in the specified format as described in README',
                          required=True)
    parser_t.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                       'characterization step (Default = training)',
                          default="training")
    parser_t.add_argument('-o', '--output', help='Output location and prefix for simulated reads (Default = simulated)',
                          default="simulated")
    parser_t.add_argument('-n', '--number', help='Number of reads to be simulated (Default = 20000)', type=int,
                          default=20000)
    parser_t.add_argument('-max', '--max_len', help='The maximum length for simulated unaligned reads. Note that this is not used for simulating aligned reads. (Default = Infinity)',
                          type=int, default=float("inf"))
    parser_t.add_argument('-min', '--min_len', help='The minimum length for simulated unaligned reads. Note that this is not used for simulating aligned reads.  (Default = 50)',
                          type=int, default=50)
    parser_t.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_t.add_argument('-hp', '--homopolymer', help='Simulate homopolymer lengths (Default = False)',
                          action='store_true', default=False)
    parser_t.add_argument('-k', '--KmerBias', help='Minimum homopolymer length to simulate homopolymer contraction and '
                                                   'expansion events in, a typical k is 6',
                          type=int, default=None)
    parser_t.add_argument('-b', '--basecaller', help='Simulate polyA tails with respect to chosen basecaller: albacore '
                                                     'or guppy',
                          choices=["albacore", "guppy"], default=None)
    parser_t.add_argument('-s', '--strandness', help='Proportion of sense sequences. Overrides the value '
                                                      'profiled in characterization stage. Should be between 0 and 1',
                          type=float, default=None)
    parser_t.add_argument('--no_model_ir', help='Ignore simulating intron retention events', action='store_false',
                          default=True)
    parser_t.add_argument('--perfect', help='Ignore profiles and simulate perfect reads', action='store_true',
                          default=False)
    parser_t.add_argument('--polya', help='Simulate polyA tails for given list of transcripts', default=None)
    parser_t.add_argument('--fastq', help='Output fastq files instead of fasta files', action='store_true',
                          default=False)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                          default=1)
    parser_t.add_argument('--uracil', help='Converts the thymine (T) bases to uracil (U) in the output fasta format',
                          action='store_true', default=False)

    parser_mg = subparsers.add_parser('metagenome', help="Run the simulator on metagenome mode")
    parser_mg.add_argument('-gl', '--genome_list', help="Reference metagenome list, tsv file, the first column is "
                                                        "species/strain name, the second column is the reference "
                                                        "genome fasta/fastq file directory", required=True)
    parser_mg.add_argument('-a', '--abun', help="Abundance list, tsv file with header, the abundance of all species in "
                                                "each sample need to sum up to 100. See example in README and provided "
                                                "config files", required=True)
    parser_mg.add_argument('-dl', '--dna_type_list',
                           help="DNA type list, tsv file, the first column is species/strain, "
                                "the second column is the chromosome name, the third column is "
                                "the DNA type: circular OR linear")
    parser_mg.add_argument('-c', '--model_prefix', help='Location and prefix of error profiles generated from '
                                                        'characterization step (Default = training)',
                           default="training")
    parser_mg.add_argument('-o', '--output',
                           help='Output location and prefix for simulated reads (Default = simulated)',
                           default="simulated")
    parser_mg.add_argument('-max', '--max_len', help='The maximum length for simulated reads (Default = Infinity)',
                           type=int, default=float("inf"))
    parser_mg.add_argument('-min', '--min_len', help='The minimum length for simulated reads (Default = 50)',
                           type=int, default=50)
    parser_mg.add_argument('-med', '--median_len', help='The median read length (Default = None), Note: this simulation'
                                                        ' is not compatible with chimeric reads simulation',
                           type=int, default=None)
    parser_mg.add_argument('-sd', '--sd_len',
                           help='The standard deviation of read length in log scale (Default = None), Note: this '
                                'simulation is not compatible with chimeric reads simulation',
                           type=float, default=None)
    parser_mg.add_argument('--seed', help='Manually seeds the pseudo-random number generator', type=int, default=None)
    parser_mg.add_argument('-hp', '--homopolymer', help=argparse.SUPPRESS,
                          action='store_true', default=False)
    parser_mg.add_argument('-k', '--KmerBias', help=argparse.SUPPRESS,
                          type=int, default=None)
    parser_mg.add_argument('-s', '--strandness', help='Percentage of antisense sequences. Overrides the value profiled '
                                                      'in characterization stage. Should be between 0 and 1',
                           type=float, default=None)
    parser_mg.add_argument('--perfect', help='Ignore error profiles and simulate perfect reads', action='store_true',
                           default=False)
    parser_mg.add_argument('--abun_var', help='Simulate random variation in abundance values, takes in two values, '
                                              'format: relative_var_low, relative_var_high, Example: -0.5 0.5)',
                           nargs='+', type=float, default=None)
    parser_mg.add_argument('--fastq', help='Output fastq files instead of fasta files', action='store_true',
                           default=False)
    parser_mg.add_argument('--chimeric', help='Simulate chimeric reads', action='store_true', default=False)
    parser_mg.add_argument('-t', '--num_threads', help='Number of threads for simulation (Default = 1)', type=int,
                           default=1)
    
    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if args.mode == "genome":
        ref_g = args.ref_g
        model_prefix = args.model_prefix
        out = args.output
        number = [args.number]
        max_len = args.max_len
        min_len = args.min_len
        median_len = args.median_len
        sd_len = args.sd_len
        chimeric = args.chimeric
        if args.seed:
            random.seed(int(args.seed))
            np.random.seed(int(args.seed))
        perfect = args.perfect
        homopolymer = args.homopolymer
        kmer_bias = args.KmerBias
        strandness = args.strandness
        dna_type = args.dna_type
        num_threads = max(args.num_threads, 1)
        fastq = args.fastq

        if homopolymer and (kmer_bias is None or kmer_bias < 0):
            print("\nPlease input proper kmer bias value >= 0 to simulate homopolymer contraction and expansion "
                  "events from\n")
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

        if median_len and sd_len and chimeric:
            sys.stderr.write("\nLognormal distributed reads cannot be chimeric!\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if max_len < min_len:
            sys.stderr.write("\nMaximum read length must be longer than Minimum read length!\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if perfect and chimeric:
            print("\nPerfect reads cannot be chimeric\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("ref_g", ref_g)
        print("model_prefix", model_prefix)
        print("out", out)
        print("number", number)
        print("perfect", perfect)
        print("homopolymer", homopolymer)
        if homopolymer:
            print("kmer_bias", kmer_bias)
        print("dna_type", dna_type)
        print("strandness", strandness)
        print("sd_len", sd_len)
        print("median_len", median_len)
        print("max_len", max_len)
        print("min_len", min_len)
        print("fastq", fastq)
        print("chimeric", chimeric)
        print("num_threads", num_threads)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        read_profile(ref_g, number, model_prefix, perfect, args.mode, strandness, dna_type=dna_type, chimeric=chimeric,
                     homopolymer=homopolymer, fastq=fastq)

        if median_len and sd_len:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Simulating read length with log-normal distribution\n")
            sys.stdout.flush()

        number_aligned = number_aligned_l[0]
        number_unaligned = number_unaligned_l[0]
        max_len = min(max_len, max_chrom)
        simulation(args.mode, out, dna_type, perfect, kmer_bias, None, max_len, min_len, num_threads, fastq, median_len,
                   sd_len, chimeric=chimeric)

    elif args.mode == "transcriptome":
        ref_g = args.ref_g
        ref_t = args.ref_t
        exp = args.exp
        model_prefix = args.model_prefix
        out = args.output
        if args.seed:
            random.seed(int(args.seed))
            np.random.seed(int(args.seed))
        number = [args.number]
        max_len = args.max_len
        min_len = args.min_len
        homopolymer = args.homopolymer
        kmer_bias = args.KmerBias
        basecaller = args.basecaller
        strandness = args.strandness
        perfect = args.perfect
        model_ir = args.no_model_ir
        dna_type = "transcriptome"
        polya = args.polya
        uracil = args.uracil
        num_threads = max(args.num_threads, 1)
        fastq = args.fastq

        if homopolymer and (kmer_bias is None or kmer_bias < 0):
            print("\nPlease input proper kmer bias value >= 0 to simulate homopolymer contraction and expansion "
                  "events from\n")
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

        if polya and basecaller is None:
            print("\nPlease input basecaller to simulate polyA tails from.\n")
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
        print("homopolymer", homopolymer)
        if homopolymer:
            print("kmer_bias", kmer_bias)
        print("model_ir", model_ir)
        print("dna_type", dna_type)
        print("strandness", strandness)
        print("max_len", max_len)
        print("min_len", min_len)
        print("uracil", uracil)
        print("polya", polya)
        if polya:
            print("basecaller", basecaller)
        print("fastq", fastq)
        print("num_threads", num_threads)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        read_profile(ref_g, number, model_prefix, perfect, args.mode, strandness, ref_t=ref_t, dna_type="linear",
                     model_ir=model_ir, polya=polya, exp=exp, homopolymer=homopolymer, fastq=fastq)

        number_aligned = number_aligned_l[0]
        number_unaligned = number_unaligned_l[0]
        max_len = min(max_len, max_chrom)
        simulation(args.mode, out, dna_type, perfect, kmer_bias, basecaller, max_len, min_len, num_threads,
                   fastq, None, None, model_ir, uracil, polya)

    elif args.mode == "metagenome":
        genome_list = args.genome_list
        abun = args.abun
        dna_type_list = args.dna_type_list
        model_prefix = args.model_prefix
        out = args.output
        max_len = args.max_len
        min_len = args.min_len
        median_len = args.median_len
        sd_len = args.sd_len
        if args.seed:
            random.seed(int(args.seed))
            np.random.seed(int(args.seed))
        perfect = args.perfect
        homopolymer = args.homopolymer
        kmer_bias = args.KmerBias
        strandness = args.strandness
        abun_var = args.abun_var
        fastq = args.fastq
        chimeric = args.chimeric
        num_threads = max(args.num_threads, 1)

        if homopolymer and (kmer_bias is None or kmer_bias < 0):
            print("\nPlease input proper kmer bias value >= 0 to simulate homopolymer contraction and expansion "
                  "events from\n")
            parser_mg.print_help(sys.stderr)
            sys.exit(1)

        if strandness and (strandness < 0 or strandness > 1):
            print("\nPlease input proper strandness value between 0 and 1\n")
            parser_mg.print_help(sys.stderr)
            sys.exit(1)

        if (median_len and not sd_len) or (sd_len and not median_len):
            sys.stderr.write("\nPlease provide both mean and standard deviation of read length!\n")
            parser_mg.print_help(sys.stderr)
            sys.exit(1)

        if median_len and sd_len and chimeric:
            sys.stderr.write("\nLognormal distributed reads cannot be chimeric!\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if max_len < min_len:
            sys.stderr.write("\nMaximum read length must be longer than Minimum read length!\n")
            parser_mg.print_help(sys.stderr)
            sys.exit(1)
            
        if perfect and chimeric:
            print("\nPerfect reads cannot be chimeric\n")
            parser_mg.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("genome_list", genome_list)
        print("abun", abun)
        print("dna_type_list", dna_type_list)
        print("model_prefix", model_prefix)
        print("out", out)
        print("perfect", perfect)
        print("strandness", strandness)
        print("sd_len", sd_len)
        print("median_len", median_len)
        print("max_len", max_len)
        print("min_len", min_len)
        print("abun_var", abun_var)
        print("fastq", fastq)
        print("chimeric", chimeric)
        print("num_threads", num_threads)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ': ' + ' '.join(sys.argv) + '\n')
        sys.stdout.flush()

        dir_name = os.path.dirname(out)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        read_profile(genome_list, [], model_prefix, perfect, args.mode, strandness, dna_type=dna_type_list, abun=abun,
                     chimeric=chimeric, homopolymer=homopolymer, fastq=fastq)

        # Add abundance variation
        global dict_abun, dict_abun_inflated
        for s in range(len(multi_dict_abun)):
            sample = list(multi_dict_abun.keys())[s]
            if abun_var:
                total_len = {}
                for species in multi_dict_abun[sample]:
                    total_len[species] = sum(seq_len[species].values())
                var_low = float(abun_var[0])
                var_high = float(abun_var[1])
                dict_abun = add_abundance_var(multi_dict_abun[sample], total_len, var_low, var_high)
            else:
                dict_abun = multi_dict_abun[sample]

            # when simulating chimeric reads, the source species for succedent segments have an inflated abundance dist
            if chimeric:
                dict_abun_inflated = {}
                for species in dict_abun:
                    dict_abun_inflated[species] = inflate_abun(dict_abun, species)

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Simulating sample " + sample + '\n')
            sys.stdout.flush()
            if median_len and sd_len:
                sys.stdout.write(
                    strftime("%Y-%m-%d %H:%M:%S") + ": Simulating read length from log-normal distribution\n")
                sys.stdout.flush()

            number_aligned = number_aligned_l[s]
            number_unaligned = number_unaligned_l[s]
            max_len = min(max_len, max(max_chrom.values()))
            simulation(args.mode, out + "_" + sample, "metagenome", perfect, kmer_bias, None, max_len,
                       min_len, num_threads, fastq, median_len, sd_len, chimeric=chimeric)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
    sys.stdout.close()


if __name__ == "__main__":
    main()

