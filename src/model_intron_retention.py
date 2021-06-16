#!/usr/bin/env python


from __future__ import with_statement
import sys
from time import strftime
import HTSeq
import pysam
from file_handler import gzopen as open

# TODO change HTSeq to Pysam

stranded = "no"


def invert_strand(iv):
    iv2 = iv.copy()
    if iv2.strand == "+":
        iv2.strand = "-"
    elif iv2.strand == "-":
        iv2.strand = "+"
    else:
        raise ValueError("Illegal strand")
    return iv2


def parse_cigar(cigar_obj):
    list_coords = []
    for co in cigar_obj:
        if co.type in ('M', '=', 'X', 'D') and co.size > 0:
            list_coords.append([co.ref_iv.chrom, co.ref_iv.start, co.ref_iv.end, co.ref_iv.strand])
    return list_coords


def intron_retention(outfile, gff_file, g_alnm, t_alnm):
    # Read intron information from GFF file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Reading intron coordinates from GFF file\n")
    gff_features = HTSeq.GFF_Reader(gff_file, end_included=True)
    features = HTSeq.GenomicArrayOfSets("auto", stranded=False)
    dict_intron_info = {}
    for feature in gff_features:
        if "transcript_id" in feature.attr:
            feature_id = feature.attr['transcript_id']
        elif "Parent" in feature.attr:  # no "if feature.type == intron" to also consider trxs without intron
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
        if feature_id not in dict_intron_info:
            dict_intron_info[feature_id] = []

        # remove "chr" from chromosome names to be constant
        if "chr" in feature.iv.chrom:
            feature.iv.chrom = feature.iv.chrom.strip("chr")

        if feature.type == "intron":
            features[feature.iv] += feature_id
            dict_intron_info[feature_id].append((feature.iv.start, feature.iv.end, feature.iv.length))

    # read primary genome alignment for each read
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read primary genome alignment for each read\n")
    dict_g_alnm = {}
    if g_alnm.lower().endswith('.bam'):
        g_alignments = pysam.AlignmentFile(g_alnm, 'rb')
    else:
        g_alignments = pysam.AlignmentFile(g_alnm)
    for r in g_alignments.fetch(until_eof=True):
        if not r.is_unmapped:
            chrom = r.reference_name
            strand = '-' if r.is_reverse else '+'
            list_coords = []
            for start, end in r.get_blocks():
                list_coords.append([chrom, start, end, strand])
            dict_g_alnm[r.query_name] = list_coords

    # read primary transcriptome alignment for each read
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read primary transcriptome alignment for each read\n")
    dict_t_alnm = {}
    if t_alnm.lower().endswith('.bam'):
        t_alignments = pysam.AlignmentFile(t_alnm, 'rb')
    else:
        t_alignments = pysam.AlignmentFile(t_alnm)
    for r in t_alignments.fetch(until_eof=True):
        if not r.is_unmapped:
            tname = r.reference_name
            if tname.startswith('ENST'):
                # an Ensembl ID
                tname = tname.split('.')[0]
            dict_t_alnm[r.query_name] = tname

    # Count the length of Intron retention events
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Calculating probabilities for each intron retention event\n")
    dict_first_intron_state = {False: 0, True: 0}
    dict_states = {(False, False): 0, (False, True): 0, (True, False): 0, (True, True): 0}
    dict_ir_info = {}
    for qname in dict_g_alnm:
        iv_seq = dict_g_alnm[qname]
        if qname in dict_t_alnm:
            primary_trx = dict_t_alnm[qname]
            if primary_trx not in dict_ir_info:
                dict_ir_info[primary_trx] = []
            list_IR_positions = []
            pos = []
            ir_info = False
            try:
                length_IR = 0
                for item in iv_seq:
                    iv = HTSeq.GenomicInterval(item[0], item[1], item[2], item[3])
                    if "chr" in iv.chrom:
                        iv.chrom = iv.chrom.strip("chr")
                    for iv2, fs2 in features[iv].steps():
                        if fs2.intersection(set([primary_trx])):
                            length_IR += iv2.length
                            pos.append(iv2.start)
                            pos.append(iv2.end)
                        else:
                            if length_IR != 0:
                                for intron in dict_intron_info[primary_trx]:
                                    if length_IR == intron[2]:
                                        list_IR_positions.append(min(pos))
                                        list_IR_positions.append(max(pos))
                                        ir_info = True
                                length_IR = 0
                                pos = []
            # TODO ??
            except UnknownChrom:
                ir_info = False
                pass

            if not ir_info:
                if primary_trx in dict_intron_info:
                    if len(dict_intron_info[primary_trx]) >= 1:  # if there is an intron
                        dict_first_intron_state[False] += 1
                        for i in range(1, len(dict_intron_info[primary_trx])):
                            dict_states[(False, False)] += 1
            else:
                # Now, go over all introns and check with the IR events
                # First we need to determine the state of first intron:
                first_intron = dict_intron_info[primary_trx][0]
                first_intron_spos = first_intron[0]
                first_intron_epos = first_intron[1]
                flag = False
                for IR_pos in list_IR_positions:
                    if first_intron_spos <= IR_pos <= first_intron_epos:
                        flag = True
                        break
                if flag:
                    dict_ir_info[primary_trx].append((first_intron_spos, first_intron_epos))
                    dict_first_intron_state[True] += 1
                    previous_state = True
                else:
                    dict_first_intron_state[False] += 1
                    previous_state = False

                # Then we will go over other introns:
                for i in range(1, len(dict_intron_info[primary_trx])):
                    intron = dict_intron_info[primary_trx][i]
                    current_state = False
                    intron_spos = intron[0]
                    intron_epos = intron[1]
                    for IR_pos in list_IR_positions:
                        if intron_spos <= IR_pos <= intron_epos:
                            current_state = True
                            dict_ir_info[primary_trx].append((intron_spos, intron_epos))
                            break
                    # print(intron_spos, intron_epos, previous_state, current_state)
                    dict_states[(previous_state, current_state)] += 1
                    previous_state = current_state

    del dict_g_alnm
    del dict_t_alnm
    # print (dict_first_intron_state)
    # print (dict_states)
    sum_first_introns = dict_first_intron_state[True] + dict_first_intron_state[False]
    sum_for_noIR = dict_states[(False, False)] + dict_states[(False, True)]
    sum_for_IR = dict_states[(True, False)] + dict_states[(True, True)]

    fout = open(outfile + "_IR_markov_model", 'w')
    fout.write("succedent\tno_IR\tIR\n")

    if sum_first_introns != 0:
        fout.write("start\t" + str(round(dict_first_intron_state[False] / float(sum_first_introns), 4)) + "\t" +
                   str(round(dict_first_intron_state[True] / float(sum_first_introns), 4)) + "\n")
    else:
        fout.write("start\t0.0\t0.0\n")

    if sum_for_noIR != 0:
        fout.write("no_IR\t" + str(round(dict_states[(False, False)] / float(sum_for_noIR), 4)) + "\t" +
                   str(round(dict_states[(False, True)] / float(sum_for_noIR), 4)) + "\n")
    else:
        fout.write("no_IR\t0.0\t0.0\n")

    if sum_for_IR != 0:
        fout.write("IR\t" + str(round(dict_states[(True, False)] / float(sum_for_IR), 4)) + "\t" +
                   str(round(dict_states[(True, True)] / float(sum_for_IR), 4)) + "\n")
    else:
        fout.write("IR\t0.0\t0.0\n")

    # output intron coordinates and information to the user:
    out_ir_info = open(outfile + "_IR_info", 'w')
    out_ir_info.write("trx_name\tintron_spos\tintron_epos\n")

    for trx in dict_ir_info:
        if len(dict_ir_info[trx]) != 0:
            lst_sorted = sorted(set(dict_ir_info[trx]))
            fstr_spos = ",".join([str(item[0]) for item in lst_sorted])
            fstr_epos = ",".join([str(item[1]) for item in lst_sorted])
            out_ir_info.write(trx + "\t" + fstr_spos + "\t" + fstr_epos + "\n")

    fout.close()
    out_ir_info.close()
