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


def primary_and_unaligned(sam_alnm_file, prefix, metagenome_list=None):
    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=True)
    if metagenome_list:
        quant_dic = {}
        for info in in_sam_file.header['SQ']:
            species_chrom = info['SN']
            species = '_'.join(species_chrom.split('_')[:-1])
            if species not in quant_dic:
                quant_dic[species] = 0  # aligned bases

    unaligned_len = []
    pos_strand = 0
    num_aligned = 0

    for aln in in_sam_file.fetch(until_eof=True):
        if not aln.is_unmapped and not aln.is_secondary and not aln.is_supplementary:
            num_aligned += 1
            out_sam_file.write(aln)
            if aln.flag == 0:
                pos_strand += 1
            if metagenome_list:
                species = '_'.join(aln.reference_name.split('_')[:-1])
                quant_dic[species] += aln.query_alignment_length
        elif aln.is_unmapped:
            unaligned_len.append(aln.query_length)

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    # Write quantification information
    if metagenome_list:
        quantification_file = open(prefix + "_quantification.tsv", 'w')
        quantification_file.write("Species\tAligned bases\tAbundance\n")
        total_bases = 0
        variation_flag = False
        for k, v in quant_dic.items():
            total_bases += v
        for k, v in quant_dic.items():
            metagenome_list[k]["real"] = v * 100 / total_bases
            quantification_file.write(k + '\t' + str(v) + '\t' + str(metagenome_list[k]["real"]) + '\n')
            if "expected" in metagenome_list[k]:
                metagenome_list[k]["variation"] = (metagenome_list[k]["real"] - metagenome_list[k]["expected"]) \
                                                  / metagenome_list[k]["expected"]
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

    return unaligned_len, strandness


def primary_and_unaligned_chimeric(sam_alnm_file, prefix, metagenome_list=None):
    """
    Function to extract alignments of reads at extremities of circular genome references
    :param sam_alnm_file: path of input SAM file
    :param prefix: prefix of output SAM file
    :param metagenome_list: metagenome dictionary containing expected abundance level in 100 scale
    :outputs: a sam file that contains: 1) primary alignments for all reads that can be aligned and 2)
              split alignments for reads that are aligned to both ends of a circular genome,
              A pickle KDE file of the gap sizes between alignment fragments
              A txt file of the Markov Model for split alignments
    :returns: an unaligned_len list, and strandness information
    """

    in_sam_file = pysam.AlignmentFile(sam_alnm_file, 'r')
    out_sam_file = pysam.AlignmentFile(prefix + "_primary.sam", 'w', template=in_sam_file, add_sam_header=True)
    chimeric_file = open(prefix + "_chimeric_info", 'w')
    gap_length = []
    chimeric_species_count = {}
    unaligned_len = []
    pos_strand = 0
    num_aligned = 0
    if metagenome_list:
        quant_dic = {}

    # extract the lengths of all reference sequences
    ref_lengths = {}
    for info in in_sam_file.header['SQ']:
        species_chrom = info['SN']
        ref_lengths[species_chrom] = info['LN']
        species = '_'.join(species_chrom.split('_')[:-1])
        chimeric_species_count[species] = [0, 0]  # the succedent segment source [same species, other species]
        if metagenome_list:
            if species not in quant_dic:
                quant_dic[species] = 0  # aligned bases

    aln_queue = []
    # parse alignment file, assuming that all alignments of a read are grouped together
    for aln in in_sam_file.fetch(until_eof=True):
        if aln.is_unmapped:
            unaligned_len.append(aln.query_length)

        elif not aln.is_secondary and not aln.is_supplementary:
            # this is a primary alignment
            num_aligned += 1

            # Define a list to store chimeric reads, each item is an interval (query_start, query_end)
            primary_direction = '+' if not aln.is_reverse else '-'
            NM_tag = int(aln.get_tag('NM'))
            primary_qstart = aln.query_alignment_start

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
                    ref_name, ref_start, direction, cigar, _ , NM_tag = supp_aln.split(',')
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
                        if len(seg["query"]) == 1 and seg["query"][0] != primary_qstart:
                            out_sam_file.write(aln)
                            # Quantification
                            if metagenome_list:
                                species = '_'.join(aln.reference_name.split('_')[:-1])
                                quant_dic[species] += aln.query_alignment_length
                            if primary_direction == '+':
                                pos_strand += 1
                            break
                        idx = [i[0] for i in sorted(enumerate(seg["query"]), key=lambda x:x[1])]
                        seg["query"].sort()
                        seg["ref"] = [seg["ref"][x] for x in idx]
                        seg["rname"] = [seg["rname"][x] for x in idx]

                        dir_added = False
                        pre_is_edge = [False, False]
                        aln_queue = [None] * len(seg["query"])
                        supplementary_to_be_added = [None] * len(seg["query"])
                        pre_species = ""
                        for i in range(len(seg["query"])):
                            interval = seg["query"][i]
                            ref_interval = seg["ref"][i]
                            is_edge = edge_checker(ref_interval[0], ref_interval[1], ref_lengths[seg["rname"][i]])
                            species = '_'.join(seg["rname"][i].split('_')[:-1])
                            if metagenome_list:
                                quant_dic[species] += seg["query"][i][1] - seg["query"][i][0]
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
                                aln_queue[i] = aln
                            else:
                                supplementary_to_be_added[i] = (seg['rname'][i], seg['query'][i][0], seg['query'][i][1],
                                                                 seg['ref'][i][0])
                            pre_is_edge = is_edge
                            pre_species = species

                        if not dir_added:
                            if seg["direction"][0] == '+':
                                pos_strand += 1
                        break

            except KeyError:
                out_sam_file.write(aln)
                # Quantification
                if metagenome_list:
                    species = '_'.join(aln.reference_name.split('_')[:-1])
                    quant_dic[species] += aln.query_alignment_length
                if primary_direction == '+':
                    pos_strand += 1

        elif aln.is_supplementary:
            qstart, qend, _, _ = cigar_parser(aln.cigarstring)
            for i in range(len(supplementary_to_be_added)):
                if (aln.reference_name, qstart, qend, aln.reference_start) == supplementary_to_be_added[i]:
                    aln_queue[i] = aln

    # output the alignments from the last read:
    for pre_aln in aln_queue:
        out_sam_file.write(pre_aln)

    in_sam_file.close()
    out_sam_file.close()

    unaligned_len = numpy.array(unaligned_len)
    strandness = float(pos_strand) / num_aligned

    # Compute the KDE of gaps
    gap_length = numpy.array(gap_length)
    gap_log = numpy.log10(gap_length + 1)
    gap_log_2d = gap_log[:, numpy.newaxis]
    kde_gap = KernelDensity(bandwidth=0.01).fit(gap_log_2d)
    joblib.dump(kde_gap, prefix + '_gap_length.pkl')

    # Write quantification information
    if metagenome_list:
        total_bases = 0
        variation_flag = False
        quantification_file = open(prefix + "_quantification.tsv", 'w')
        quantification_file.write("Species\tAligned bases\tAbundance\n")

        for k, v in quant_dic.items():
            total_bases += v
        for k, v in quant_dic.items():
            metagenome_list[k]["real"] = v * 100 / total_bases
            quantification_file.write(k + '\t' + str(v) + '\t' + str(metagenome_list[k]["real"]) + '\n')
            if "expected" in metagenome_list[k]:
                metagenome_list[k]["variation"] = (metagenome_list[k]["real"] - metagenome_list[k]["expected"])\
                                                  / metagenome_list[k]["expected"]
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

    # Compute chimeric information
    mean_segments = (len(gap_length) + num_aligned) / num_aligned
    chimeric_file.write("Mean segments for each aligned read:\t" + str(mean_segments) + '\n')
    if metagenome_list:
        chimeric_file.write("Shrinkage rate (beta):\t" + str(median(beta_list)))
    chimeric_file.close()

    return unaligned_len, strandness

