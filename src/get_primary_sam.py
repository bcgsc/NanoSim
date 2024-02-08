#!/usr/bin/env python

from __future__ import with_statement
import pysam
import sys
import numpy
import re
import joblib
from sklearn.neighbors import KernelDensity
from head_align_tail_dist import edge_checker
from statistics import median
from time import strftime
from file_handler import gzopen as open


def cigar_parser(cigar):
    match = re.findall(r'(\d+)(\w)', cigar)
    qstart = int(match[0][0]) if match[0][1] in ('S', 'H') else 0
    qlen = 0
    rlen = 0

    for item in match:
        if item[1] == 'M':
            qlen += int(item[0])
            rlen += int(item[0])
        elif item[1] == 'I':
            qlen += int(item[0])
        elif item[1] == 'D':
            rlen += int(item[0])
    qend = qstart + qlen
    return qstart, qend, qlen, rlen


def not_overlap(interval, interval_lst, interval_name=None, interval_name_list=None, overlap_base=10):
    # interval: (start, end)
    # overlap_base: a relaxed threshold reserved for a small stretch of bases shared by both alignments
    for i in range(len(interval_lst)):
        if interval[0] < interval_lst[i][1] - overlap_base and interval[1] - overlap_base > interval_lst[i][0]:
            if interval_name is None or (interval_name is not None and interval_name == interval_name_list[i]):
                return False
    return True


def EM_meta(read_list, all_species):
    """
    :param read_list: {(read1, interval1): [species1, species2], (read1, interval2): [species1, species2] ...}
    :return: abundance_list: {species1: abundance1, species2: abundance2, ...}
    """
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Starting EM for quantification\n")
    base_count_unique = dict.fromkeys(all_species, 0)
    multi_mapping_list = {}
    total_base = 0
    for read, species_list in read_list.items():
        length = read[1][1] - read[1][0]
        total_base += length
        if len(species_list) == 1:
            base_count_unique[species_list[0]] += length
        else:
            multi_mapping_list[(read[0], length)] = species_list
    abundance_list = {species:  100 / len(all_species) for species in all_species}

    diff = 100 * len(all_species)
    for iter in range(0, 100):
        # Expectation
        base_count_tmp = {k: v for k, v in base_count_unique.items()}
        for read, species_list in multi_mapping_list.items():
            length = read[1]
            total_abun_species = sum([abundance_list[species] for species in species_list])
            perc_for_species = [abundance_list[species] / total_abun_species for species in species_list]
            for i in range(len(species_list)):
                base_count_tmp[species_list[i]] += length * perc_for_species[i]

        # Maximization
        abundance_list_tmp = {species:  bases * 100 / total_base for species, bases in base_count_tmp.items()}
        diff_tmp = sum([abs(abundance_list_tmp[i] - abundance_list[i]) for i in abundance_list.keys()])
        if iter % 10 == 0:
            sys.stdout.write("Iteration: " + str(iter) + ", Diff: " + str(diff_tmp) + '\n')
            sys.stdout.flush()

        abundance_list = abundance_list_tmp
        diff_thres = min(abundance_list.values()) * 0.01
        if diff_tmp <= diff_thres or diff - diff_tmp < diff_thres:
            break
        diff = diff_tmp

    return abundance_list


def EM_trans(read_list, all_trans, normalize):
    """
    :param read_list: {(read1, interval1): [trans1, trans2]), (read1, interval2), [trans1, trans2]) ...}
    :return: tpm_list: {species1: tpm1, species2: tpm2, ...}
    """
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Starting EM for quantification\n")
    read_count_unique = dict.fromkeys(all_trans, 0)
    multi_mapping_list = {}
    total_reads = 0
    for read, trans_list in read_list.items():
        total_reads += 1
        if len(trans_list) == 1:
            read_count_unique[trans_list[0]] += 1
        else:
            multi_mapping_list[read[0]] = trans_list
    abundance_list = {trans: 100 / len(all_trans) for trans in all_trans}

    diff = 100 * len(all_trans)
    for iter in range(0, 1000):
        # Expectation
        read_count_tmp = {k: v for k, v in read_count_unique.items()}
        for read, trans_list in multi_mapping_list.items():
            all_trans_in_read = sum([abundance_list[trans] for trans in trans_list])
            perc_for_trans = [abundance_list[trans] / all_trans_in_read for trans in trans_list]
            for i in range(len(trans_list)):
                read_count_tmp[trans_list[i]] += perc_for_trans[i]

        # Maximization
        abundance_list_tmp = {trans: reads * 100 / total_reads for trans, reads in read_count_tmp.items()}
        diff_tmp = sum([abs(abundance_list_tmp[i] - abundance_list[i]) for i in abundance_list.keys()])
        if iter % 50 == 0:
            sys.stdout.write("Iteration: " + str(iter) + ", Diff: " + str(diff_tmp) + '\n')
            sys.stdout.flush()

        abundance_list = abundance_list_tmp
        diff_thres = min(abundance_list.values()) * 0.001
        if diff_tmp <= diff_thres or diff - diff_tmp < diff_thres:
            break
        diff = diff_tmp

    # calculate TPM
    tpm_list = {}
    total_rpk = 0
    if normalize:
        for trans, count in read_count_tmp.items():
            total_rpk += count / all_trans[trans] * 1e3  # read per kilobase
    else:
        total_rpk = sum(read_count_tmp.values())
    for trans, count in read_count_tmp.items():
        rpk = count / all_trans[trans] * 1e3 if normalize else count
        tpm = rpk * 1e6 / total_rpk
        tpm_list[trans] = (count, tpm)

    return tpm_list


def primary_and_unaligned(sam_alnm_file, prefix, metagenome_list=None, fastq=False):
    in_sam_file = pysam.AlignmentFile(sam_alnm_file)
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.bam", 'wb', template=in_sam_file, add_sam_header=True)
    if metagenome_list:
        quant_dic = {}

    unaligned_len = []
    pos_strand = 0
    num_aligned = 0

    unaligned_bq = []

    all_species = {}
    for info in in_sam_file.header['SQ']:
        species_chrom = info['SN']
        species = '_'.join(species_chrom.split('_')[:-1])
        all_species[species] = 0

    for aln in in_sam_file.fetch(until_eof=True):
        if not aln.is_unmapped and not aln.is_secondary and not aln.is_supplementary:
            num_aligned += 1
            out_sam_file.write(aln)
            if aln.flag == 0:
                pos_strand += 1
            if metagenome_list:
                species = '_'.join(aln.reference_name.split('_')[:-1])
                quant_dic[(aln.query_name, (aln.query_alignment_start, aln.query_alignment_end))] = [species]
        elif aln.is_unmapped:
            unaligned_len.append(aln.query_length)
            if fastq and aln.query_alignment_qualities:
                unaligned_bq += aln.query_alignment_qualities.tolist()
        elif aln.is_secondary or aln.is_supplementary:
            if metagenome_list:
                qstart, qend, _, _ = cigar_parser(aln.cigarstring)
                if (aln.query_name, (qstart, qend)) in quant_dic:
                    species = '_'.join(aln.reference_name.split('_')[:-1])
                    quant_dic[(aln.query_name, (qstart, qend))].append(species)

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    # Write quantification information
    if metagenome_list:
        variation_flag = False
        quantification_file = open(prefix + "_quantification.tsv", 'w')
        quantification_file.write("Species\tAbundance\n")

        abundance_list = EM_meta(quant_dic, all_species)
        for k, v in abundance_list.items():
            quantification_file.write(k + '\t' + str(v) + '\n')
            metagenome_list[k]["real"] = v
            if "expected" in metagenome_list[k]:
                metagenome_list[k]["variation"] = (v - metagenome_list[k]["expected"]) / metagenome_list[k]["expected"]
                variation_flag = True

        quantification_file.close()

        if variation_flag:
            variations = [v["variation"] for v in metagenome_list.values()]
            var_low = min(variations)
            var_high = max(variations)
            var_median = median([abs(v) for v in variations])
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Variation in abundance (low, high, abs(median)): "
                                                             "" + str(var_low) + ',' + str(var_high) + ',' + \
                             str(var_median) + '\n')
            sys.stdout.flush()

        quantification_file.close()

    return unaligned_len, strandness, unaligned_bq


def primary_and_unaligned_chimeric(sam_alnm_file, prefix, metagenome_list=None, q_mode=False, normalize=True, fastq=False):
    """
    Function to extract alignments of reads at extremities of circular genome references
    :param sam_alnm_file: path of input SAM file
    :param prefix: prefix of output BAM file
    :param metagenome_list: metagenome dictionary containing expected abundance level in 100 scale
    :param q_mode: whether it is quantification only mode
    :outputs: a sam file that contains: 1) primary alignments for all reads that can be aligned and 2)
              split alignments for reads that are aligned to both ends of a circular genome,
              A pickle KDE file of the gap sizes between alignment fragments
              A txt file of the Markov Model for split alignments
    :returns: an unaligned_len list, and strandness information
    """

    is_trans = False

    in_sam_file = pysam.AlignmentFile(sam_alnm_file)
    if metagenome_list:
        quant_dic = {}
        is_trans = True if 'tpm' in metagenome_list else False

    if not q_mode:
        out_sam_file = pysam.AlignmentFile(prefix + "_primary.bam", 'wb', template=in_sam_file, add_sam_header=True)
        chimeric_file = open(prefix + "_chimeric_info", 'w')
    gap_length = []
    chimeric_species_count = {}
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
    all_species = {}
    unaligned_bq = []

    # extract the lengths of all reference sequences
    ref_lengths = {}
    for info in in_sam_file.header['SQ']:
        species_chrom = info['SN']
        ref_lengths[species_chrom] = info['LN']
        if metagenome_list and is_trans:
            species = species_chrom
        else:
            species = '_'.join(species_chrom.split('_')[:-1])
        all_species[species] = info['LN']
        chimeric_species_count[species] = [0, 0]  # the succedent segment source [same species, other species]

    if not q_mode:
        aln_queue = []
    # parse alignment file, assuming that all alignments of a read are grouped together
    for aln in in_sam_file.fetch(until_eof=True):
        if aln.is_unmapped:
            unaligned_len.append(aln.query_length)
            if fastq and aln.query_alignment_qualities:
                unaligned_bq += aln.query_alignment_qualities.tolist()

        elif not aln.is_secondary and not aln.is_supplementary:
            # this is a primary alignment
            num_aligned += 1

            # Define a list to store chimeric reads, each item is an interval (query_start, query_end)
            primary_direction = '+' if not aln.is_reverse else '-'
            NM_tag = int(aln.get_tag('NM'))
            primary_qstart = aln.query_alignment_start

            if not q_mode:
                supplementary_to_be_added = []
                # output the alignments from the previous read first:
                for pre_aln in aln_queue:
                    out_sam_file.write(pre_aln)
                aln_queue = []

            try:
                supplementary_aln_list = aln.get_tag('SA').split(';')
                compatible_list = [{"query": [(aln.query_alignment_start, aln.query_alignment_end)],
                                    "ref": [(aln.reference_start, aln.reference_end)],
                                    "score": aln.query_alignment_length - NM_tag,
                                    "rname": [aln.reference_name],
                                    "direction": [primary_direction]}]
                # parse supplementary alignments and add to compatible lists
                for supp_aln in supplementary_aln_list[:-1]:
                    ref_name, ref_start, direction, cigar, _, NM_tag = supp_aln.split(',')
                    ref_start = int(ref_start) - 1
                    NM_tag = int(NM_tag)
                    qstart, qend, qlen, rlen = cigar_parser(cigar)

                    added = False
                    for seg in compatible_list:
                        if not_overlap((qstart, qend), seg["query"]) and \
                                not_overlap((ref_start, ref_start + rlen), seg["ref"], ref_name, seg['rname']):
                            seg["query"].append((qstart, qend))
                            seg["ref"].append((ref_start, ref_start + rlen))
                            seg["score"] += (qlen - NM_tag)
                            seg["rname"].append(ref_name)
                            seg["direction"].append(direction)
                            added = True

                    if not added:
                        compatible_list.append({"query": [(qstart, qend)],
                                                "ref": [(ref_start, ref_start + rlen)],
                                                "score": qlen - NM_tag,
                                                "rname": [ref_name],
                                                "direction": [direction]})

                max_score = max([x["score"] for x in compatible_list])
                for seg in compatible_list:
                    if seg["score"] == max_score:
                        # if the only alignment picked is not primary alignment, change it to primary alignment
                        if len(seg["query"]) == 1 and seg["query"][0][0] != primary_qstart:
                            if not q_mode:
                                out_sam_file.write(aln)
                            # Quantification
                            if metagenome_list and is_trans:
                                species = aln.reference_name
                                quant_dic[(aln.query_name, (seg["query"][0][0], seg["query"][0][1]))] = [species]
                            elif metagenome_list:
                                species = '_'.join(aln.reference_name.split('_')[:-1])
                                quant_dic[(aln.query_name, (seg["query"][0][0], seg["query"][0][1]))] = [species]
                            if primary_direction == '+':
                                pos_strand += 1
                            break
                        idx = [i[0] for i in sorted(enumerate(seg["query"]), key=lambda x: x[1])]
                        seg["query"].sort()
                        seg["ref"] = [seg["ref"][x] for x in idx]
                        seg["rname"] = [seg["rname"][x] for x in idx]

                        dir_added = False
                        pre_is_edge = [False, False]
                        if not q_mode:
                            aln_queue = [None] * len(seg["query"])
                            supplementary_to_be_added = [None] * len(seg["query"])
                        pre_species = ""
                        for i in range(len(seg["query"])):
                            interval = seg["query"][i]
                            ref_interval = seg["ref"][i]
                            is_edge = edge_checker(ref_interval[0], ref_interval[1], ref_lengths[seg["rname"][i]])
                            if metagenome_list:
                                species = seg["rname"][i] if is_trans else '_'.join(seg["rname"][i].split('_')[:-1])
                                quant_dic[(aln.query_name, (interval[0], interval[1]))] = [species]
                            # Record the gap
                            if i > 0:
                                # The gap between split alignments of circular read is not recorded
                                if seg["rname"][i] == seg["rname"][i - 1] and \
                                        ((pre_is_edge[0] and is_edge[1]) or (pre_is_edge[1] and is_edge[0])):
                                    pass  # circular reads
                                else:
                                    gap = max(0, interval[0] - seg["query"][i - 1][1])  # Change negative gaps size to 0
                                    gap_length.append(gap)
                                    if species == pre_species:
                                        chimeric_species_count[pre_species][0] += 1
                                    else:
                                        chimeric_species_count[pre_species][1] += 1
                            if interval[0] == primary_qstart:
                                dir_added = True
                                if primary_direction == '+':
                                    pos_strand += 1
                                if not q_mode:
                                    aln_queue[i] = aln
                            else:
                                if not q_mode:
                                    supplementary_to_be_added[i] = (seg['rname'][i], seg['query'][i][0],
                                                                    seg['query'][i][1], seg['ref'][i][0])
                            pre_is_edge = is_edge
                            pre_species = species

                        if not dir_added:
                            if seg["direction"][0] == '+':
                                pos_strand += 1
                        break

            except KeyError:
                if not q_mode:
                    out_sam_file.write(aln)
                # Quantification
                if metagenome_list:
                    species = aln.reference_name if is_trans else '_'.join(aln.reference_name.split('_')[:-1])
                    quant_dic[(aln.query_name, (aln.query_alignment_start, aln.query_alignment_end))] = [species]
                if primary_direction == '+':
                    pos_strand += 1

        elif aln.is_supplementary or aln.is_secondary:
            qstart, qend, _, _ = cigar_parser(aln.cigarstring)
            if not q_mode:
                for i in range(len(supplementary_to_be_added)):
                    if (aln.reference_name, qstart, qend, aln.reference_start) == supplementary_to_be_added[i]:
                        aln_queue[i] = aln
            if metagenome_list and (aln.query_name, (qstart, qend)) in quant_dic:
                species = aln.reference_name if is_trans else '_'.join(aln.reference_name.split('_')[:-1])
                quant_dic[(aln.query_name, (qstart, qend))].append(species)

    if not q_mode:
        # output the alignments from the last read:
        for pre_aln in aln_queue:
            out_sam_file.write(pre_aln)

    in_sam_file.close()

    # Write quantification information
    if is_trans:
        quantification_file = open(prefix + "_quantification.tsv", 'w')
        tpm_list = EM_trans(quant_dic, all_species, normalize)
        quantification_file.write("ID\tcount\tTPM\n")
        for trans, info in tpm_list.items():
            quantification_file.write(trans + '\t' + str(info[0]) + '\t' + str(info[1]) + '\n')
        quantification_file.close()

    elif metagenome_list:
        quantification_file = open(prefix + "_quantification.tsv", 'w')
        variation_flag = False
        quantification_file.write("Species\tAbundance\n")

        abundance_list = EM_meta(quant_dic, all_species)
        for k, v in abundance_list.items():
            quantification_file.write(k + '\t' + str(v) + '\n')
            metagenome_list[k]["real"] = v
            if "expected" in metagenome_list[k]:
                metagenome_list[k]["variation"] = (v - metagenome_list[k]["expected"]) / metagenome_list[k]["expected"]
                variation_flag = True

        quantification_file.close()

        if variation_flag:
            variations = [v["variation"] for v in metagenome_list.values()]
            var_low = min(variations)
            var_high = max(variations)
            var_median = median([abs(v) for v in variations])
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Variation in abundance (low, high, abs(median)): "
                                                             "" + str(var_low) + ',' + str(var_high) + ',' + \
                             str(var_median) + '\n')
            sys.stdout.flush()

        # inflated_prob_for_one_species + beta * prob_for_other_species = 1 for chimeric reads
        beta_list = []
        for species, counts in chimeric_species_count.items():
            original_prob = metagenome_list[species]["real"]
            other_prob = 100 - original_prob
            if counts[0] + counts[1] == 0:
                continue
            beta_list.append(counts[1] / (counts[0] + counts[1]) * 100 / other_prob)

    strandness = float(pos_strand) / num_aligned
    if q_mode:
        return [], strandness

    out_sam_file.close()
    unaligned_len = numpy.array(unaligned_len)

    # Compute the KDE of gaps
    gap_length = numpy.array(gap_length)
    gap_log = numpy.log10(gap_length + 1)
    gap_log_2d = gap_log[:, numpy.newaxis]
    kde_gap = KernelDensity(bandwidth=0.01).fit(gap_log_2d)
    joblib.dump(kde_gap, prefix + '_gap_length.pkl')

    # Compute chimeric information
    mean_segments = (len(gap_length) + num_aligned) / num_aligned
    chimeric_file.write("Mean segments for each aligned read:\t" + str(mean_segments) + '\n')
    if metagenome_list and not is_trans:
        chimeric_file.write("Shrinkage rate (beta):\t" + str(median(beta_list)))
    chimeric_file.close()

    return unaligned_len, strandness, unaligned_bq
