#!/usr/bin/env python

"""
@author: Chen Yang & Saber HafezQorani
This script generates read profiles for Oxford Nanopore 2D reads (genomic and transcriptome).
"""


from __future__ import print_function
from __future__ import with_statement
from textwrap import dedent
from subprocess import call
from time import strftime
import sys
import os
import re
import argparse
import numpy
from sklearn.neighbors import KernelDensity
from sklearn.externals import joblib
import head_align_tail_dist as align
import get_besthit_maf
import get_primary_sam
import besthit_to_histogram as error_model
import model_fitting
import model_intron_retention as model_ir


PYTHON_VERSION = sys.version_info
VERSION = "2.5.0"
PRORAM = "NanoSim"
AUTHOR = "Chen Yang, Saber Hafezqorani (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca; shafezqorani@bcgsc.ca"


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


def align_transcriptome(in_fasta, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g, post=True):
    if t_alnm == '':
        if aligner == "minimap2":
            t_alnm = prefix + "_transcriptome_alnm.sam"
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
            call("minimap2 --cs -ax map-ont -t " + num_threads + " " + ref_t + " " + in_fasta + " > " + t_alnm,
                 shell=True)

        elif aligner == "LAST":
            t_alnm = prefix + "_transcriptome_alnm.maf"
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference transcriptome\n")
            call("lastdb ref_transcriptome " + ref_t, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_transcriptome " + in_fasta + " > " + t_alnm, shell=True)

    if g_alnm == '':
        if aligner == "minimap2":
            g_alnm = prefix + "_genome_alnm.sam"
            # Alignment to reference genome
            # [EDIT] I may change the options for minimap2 when dealing with cDNA and dRNA reads.
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference genome\n")
            call("minimap2 --cs -ax splice -t " + num_threads + " " + ref_g + " " + in_fasta + " > " + g_alnm, shell=True)

        elif aligner == "LAST":
            g_alnm = prefix + "_genome_alnm.maf"
            # Alignment to reference genome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference genome\n")
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " > " + g_alnm, shell=True)

    if not post:
        return t_alnm, g_alnm

    # post-process
    t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
    t_alnm_ext = t_alnm_ext[1:]
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing transcriptome alignment file: " + t_alnm_ext + '\n')
    if t_alnm_ext == "maf":
        processed_maf_t = prefix + "_transcriptome_alnm_processed.maf"
        call("grep '^s ' " + t_alnm + " > " + processed_maf_t, shell=True)
        unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf_t, prefix + "_transcriptome")
    elif t_alnm_ext == "sam":
        unaligned_length, strandness = get_primary_sam.primary_and_unaligned(t_alnm, prefix + "_transcriptome")

    g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
    g_alnm_ext = g_alnm_ext[1:]
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing genome alignment file: " + g_alnm_ext + '\n')
    if g_alnm_ext == "maf":
        processed_maf = prefix + "_processed.maf"
        call("grep '^s ' " + g_alnm + " > " + processed_maf, shell=True)
        get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf, prefix + "_genome")
    elif g_alnm_ext == "sam":
        get_primary_sam.primary_and_unaligned(g_alnm, prefix + "_genome")

    return t_alnm_ext, unaligned_length, g_alnm, t_alnm, strandness


def align_genome(in_fasta, prefix, aligner, num_threads, g_alnm, ref_g, post=True):
    # if an alignment file is not provided
    if g_alnm == '':
        if aligner == "minimap2":  # Align with minimap2 by default
            g_alnm = prefix + "_genome_alnm.sam"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2\n")
            call("minimap2 --cs --MD -ax map-ont -t " + num_threads + " " + ref_g + " " + in_fasta + " > " + g_alnm,
                 shell=True)

        elif aligner == "LAST":
            g_alnm = prefix + "_genome_alnm.maf"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " " + g_alnm, shell=True)

    if not post:
        return

    # post-process
    pre, file_ext = os.path.splitext(g_alnm)
    file_extension = file_ext[1:]
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing alignment file: " + file_extension + "\n")
    if file_extension == "maf":
        processed_maf = prefix + "_processed.maf"
        call("grep '^s ' " + g_alnm + " > " + processed_maf, shell=True)

        # get best hit and unaligned reads
        unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf, prefix)

    elif file_extension == "sam":
        # get the primary alignments and define unaligned reads.
        unaligned_length, strandness = get_primary_sam.primary_and_unaligned(g_alnm, prefix)

    return file_extension, unaligned_length, strandness


def add_intron(annot, prefix):
    # Read the annotation GTF/GFF3 file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
    # If gtf provided, convert to GFF3 (gt gtf_to_gff3)
    annot_filename, annot_file_extension = os.path.splitext(annot)
    annot_file_extension = annot_file_extension[1:]
    if annot_file_extension.upper() == "GTF":
        call("gt gtf_to_gff3 -tidy -force -o " + prefix + ".gff3 " + annot, shell=True)
    else:
        call("cp " + annot + ' ' + prefix + '.gff3', shell=True)

    # TODO Check whether the gff3 file already has intron infomation

    # Next, add intron info into gff3:
    call("gt gff3 -tidy -retainids -checkids -addintrons -sort -force -o " + prefix + "_added_intron_temp.gff3 " +
         prefix + ".gff3", shell=True)

    # Inherit "transcript_id" information for intron features from exon info
    script_path = os.path.realpath(__file__)
    script_dir = os.path.dirname(script_path)

    call("gt " + script_dir + "/bequeath.lua transcript_id < " + prefix + "_added_intron_temp.gff3 > " + prefix +
         "_added_intron_final.gff3", shell=True)
    call("rm " + prefix + "_added_intron_temp.gff3", shell=True)


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Read characterization step
        -----------------------------------------------------------
        Given raw ONT reads, reference genome and/or transcriptome,
        learn read features and output error profiles
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='NanoSim ' + VERSION)
    subparsers = parser.add_subparsers(dest='mode', description=dedent('''
        There are four modes in read_analysis.
        For detailed usage of each mode:
            read_analysis.py mode -h
        -------------------------------------------------------
        '''))

    parser_g = subparsers.add_parser('genome', help="Run the simulator on genome mode")
    parser_g.add_argument('-i', '--read', help='Input read for training', required=True)
    parser_g.add_argument('-rg', '--ref_g', help='Reference genome, not required if genome alignment file is provided',
                          default='')
    parser_g.add_argument('-a', '--aligner', help='The aligner to be used, minimap2 or LAST (Default = minimap2)',
                          choices=['minimap2', 'LAST'], default='minimap2')
    parser_g.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format (optional)', default='')
    parser_g.add_argument('-o', '--output', help='The location and prefix of outputting profiles (Default = training)',
                          default='training')
    parser_g.add_argument('--no_model_fit', help='Disable model fitting step', action='store_false', default=True)
    parser_g.add_argument('-t', '--num_threads', help='Number of threads for alignment and model fitting (Default = 1)',
                          type=int, default=1)

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode")
    parser_t.add_argument('-i', '--read', help='Input read for training', required=True)
    parser_t.add_argument('-rg', '--ref_g', help='Reference genome', required=True)
    parser_t.add_argument('-rt', '--ref_t', help='Reference Transcriptome', required=True)  # ?
    parser_t.add_argument('-annot', '--annotation', help='Annotation file in ensemble GTF/GFF formats, '
                                                         'required for intron retention detection', default='')
    parser_t.add_argument('-a', '--aligner', help='The aligner to be used: minimap2 or LAST (Default = minimap2)',
                          choices=['minimap2', 'LAST'], default='minimap2')
    parser_t.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format (optional)', default='')
    parser_t.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in sam or maf format (optional)',
                          default='')
    parser_t.add_argument('-o', '--output', help='The location and prefix of outputting profiles (Default = training)',
                          default='training')
    parser_t.add_argument('--no_model_fit', help='Disable model fitting step', action='store_false', default=True)
    parser_t.add_argument('--no_intron_retention', help='Disable Intron Retention analysis', action='store_false',
                          default=True)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads for alignment and model fitting (Default = 1)',
                          type=int, default=1)

    parser_e = subparsers.add_parser('quantify', help="Quantify expression profile of transcripts")
    parser_e.add_argument('-i', '--read', help='Input reads for quantification', required=True)
    parser_e.add_argument('-rt', '--ref_t', help='Reference Transcriptome', required=True)
    parser_e.add_argument('-o', '--output', help='The location and prefix of outputting profile (Default = expression)',
                          default='expression')
    parser_e.add_argument('-t', '--num_threads', help='Number of threads for alignment (Default = 1)', type=int,
                          default=1)

    parser_ir = subparsers.add_parser('detect_ir', help="Detect Intron Retention events using the alignment file")
    parser_ir.add_argument('-annot', '--annotation', help='Annotation file in ensemble GTF/GFF formats', required=True)
    parser_ir.add_argument('-i', '--read', help='Input read for training, not required if alignment files are provided',
                           default='')
    parser_ir.add_argument('-rg', '--ref_g', help='Reference genome, not required if genome alignment file is provided',
                           default='')
    parser_ir.add_argument('-rt', '--ref_t', help='Reference Transcriptome, not required if transcriptome alignment '
                                                  'file is provided', default='')
    parser_ir.add_argument('-a', '--aligner', help='The aligner to be used: minimap2 or LAST (Default = minimap2)',
                           choices=['minimap2', 'LAST'], default='minimap2')
    parser_ir.add_argument('-o', '--output', help='The output name and location', required=False, default='ir_info')
    parser_ir.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format (optional)', default='')
    parser_ir.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in sam or maf format (optional)',
                           default='')
    parser_ir.add_argument('-t', '--num_threads', help='Number of threads for alignment (Default = 1)', type=int,
                           default=1)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # parse quantify mode arguments
    if args.mode == "quantify":
        infile = args.read
        ref_t = args.ref_t
        prefix = args.output
        num_threads = str(max(args.num_threads, 1))

        print("\nrunning the code with following parameters:\n")
        print("infile", infile)
        print("ref_t", ref_t)
        print("prefix", prefix)
        print("num_threads", num_threads)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # Quantifying the transcript abundance from input read
        sys.stdout.write('Quantifying transcripts abundance: \n')
        map_file = prefix + '_mapping.paf'
        call("minimap2 -t " + str(num_threads) + " -x map-ont -p0 " + ref_t + " " + infile + " > " + map_file,
             shell=True)

        # Get the script path
        script_path = os.path.realpath(__file__)
        script_dir = os.path.dirname(script_path)
        out_file = prefix + '_abundance.tsv'
        call("python " + script_dir + "/nanopore_transcript_abundance.py -i " + map_file + " > " + out_file, shell=True)
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
        return

    # parse detect_ir mode arguments
    if args.mode == "detect_ir":
        annot = args.annotation
        infile = args.read
        prefix = args.output
        aligner = args.aligner
        ref_g = args.ref_g
        ref_t = args.ref_t
        g_alnm = args.g_alnm
        t_alnm = args.t_alnm
        num_threads = str(max(args.num_threads, 1))

        if g_alnm == '' and ref_g == '':
            print("Please supply a reference genome or genome alignment file\n")
            parser_ir.print_help(sys.stderr)
            sys.exit(1)

        if t_alnm == '' and ref_t == '':
            print("Please supply a reference transcriptome or transcriptome alignment file\n")
            parser_ir.print_help(sys.stderr)
            sys.exit(1)

        # check validity of parameters
        if g_alnm != '':
            pre, file_ext = os.path.splitext(g_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam']:
                print("Please specify an acceptable alignment format! (.maf or .sam)\n")
                parser_ir.print_help(sys.stderr)
                sys.exit(1)
        if t_alnm != '':
            pre, file_ext = os.path.splitext(t_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam']:
                print("Please specify an acceptable alignment format! (.maf or .sam)\n")
                parser_ir.print_help(sys.stderr)
                sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("annot", annot)
        print("infile", infile)
        print("aligner", aligner)
        print("ref_g", ref_g)
        print("ref_t", ref_t)
        print("g_alnm", g_alnm)
        print("t_alnm", t_alnm)
        print("prefix", prefix)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # Alignment if maf/sam file not provided, and post process them to include only primary alignments
        alnm_ext, unaligned_length, g_alnm, t_alnm, strandness \
            = align_transcriptome(infile, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g)

        # Add introns to annotation file
        add_intron(annot, prefix)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
        model_ir.intron_retention(prefix, prefix + "_added_intron_final.gff3", g_alnm, t_alnm)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
        return

    if args.mode == "genome":
        infile = args.read
        ref_g = args.ref_g
        aligner = args.aligner
        g_alnm = args.g_alnm
        prefix = args.output
        num_threads = str(max(args.num_threads, 1))
        model_fit = args.no_model_fit

        # check validity of parameters
        if g_alnm != '':
            pre, file_ext = os.path.splitext(g_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam']:
                print("Please specify an acceptable alignment format! (.maf or .sam)\n")
                parser_g.print_help(sys.stderr)
                sys.exit(1)
        if g_alnm == '' and ref_g == '':
            print("Please supply a reference genome or genome alignment file\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        print("\nRunning the code with following parameters:\n")
        print("infile", infile)
        print("ref_g", ref_g)
        print("aligner", aligner)
        print("g_alnm", g_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("model_fit", model_fit)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process\n")
        in_fasta = prefix + "_processed.fasta"
        processed_fasta = open(in_fasta, 'w')

        # Replace spaces in sequence headers with dashes to create unique header for each read
        with open(infile, 'r') as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                processed_fasta.write('>' + chr_name + '\n' + seqS + '\n')
        processed_fasta.close()

        alnm_ext, unaligned_length, strandness = align_genome(in_fasta, prefix, aligner, num_threads, g_alnm, ref_g)

        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        num_aligned = align.head_align_tail(prefix, alnm_ext, args.mode)

    if args.mode == "transcriptome":
        infile = args.read
        ref_g = args.ref_g
        ref_t = args.ref_t
        annot = args.annotation
        aligner = args.aligner
        g_alnm = args.g_alnm
        t_alnm = args.t_alnm
        prefix = args.output
        num_threads = str(max(args.num_threads, 1))
        model_fit = args.no_model_fit
        ir = args.no_intron_retention

        if ir and g_alnm == '' and ref_g == '':
            print("For intron retention function, please supply a reference genome or genome alignment file\n")
            parser_ir.print_help(sys.stderr)
            sys.exit(1)

        if t_alnm == '' and ref_t == '':
            print("Please supply a reference transcriptome or transcriptome alignment file\n")
            parser_ir.print_help(sys.stderr)
            sys.exit(1)

        if g_alnm != '' and t_alnm != '':
            g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
            t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
            g_alnm_ext = g_alnm_ext[1:]
            t_alnm_ext = t_alnm_ext[1:]
            if g_alnm_ext != t_alnm_ext or g_alnm_ext not in ['maf', 'sam']:
                print("\nPlease provide both alignments in a same format: sam OR maf\n")
                parser_t.print_help(sys.stderr)
                sys.exit(1)
            # Development: model IR using MAF alignment formats as well
            if g_alnm_ext == t_alnm_ext == "maf" and ir:
                print("\nThe intron retention only works with sam alignment files for now. Thanks\n")
                parser_t.print_help(sys.stderr)
                sys.exit(1)

        if ir and (ref_g == '' or annot == ''):
            print("\nPlease also input reference genome and annotation file for Intron Retention modeling\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        print("\nrunning the code with following parameters:\n")
        print("infile", infile)
        print("ref_g", ref_g)
        print("ref_t", ref_t)
        print("annot", annot)
        print("aligner", aligner)
        print("g_alnm", g_alnm)
        print("t_alnm", t_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("model_fit", model_fit)
        print("intron_retention", ir)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")
        in_fasta = prefix + "_processed.fasta"
        processed_fasta = open(in_fasta, 'w')
        with open(infile, 'r') as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                processed_fasta.write('>' + chr_name + '\n' + seqS + '\n')
        processed_fasta.close()

        # Read the length of reference transcripts from the reference transcriptome
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read the length of reference transcripts \n")
        dict_ref_len = {}
        with open(ref_t) as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                dict_ref_len[chr_name] = len(seqS)

        alnm_ext, unaligned_length, g_alnm, t_alnm, strandness = \
            align_transcriptome(in_fasta, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g)

        if ir:
            # Add introns to annotation file
            add_intron(annot, prefix)

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
            model_ir.intron_retention(prefix, prefix + "_added_intron_final.gff3", g_alnm, t_alnm)

        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        num_aligned = align.head_align_tail(prefix + "_transcriptome", alnm_ext, args.mode, dict_ref_len)

    # strandness of the aligned reads
    strandness_rate = open(prefix + "_strandness_rate", 'w')
    strandness_rate.write("strandness:\t" + str(round(strandness, 3)))
    strandness_rate.close()

    # Length distribution of unaligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Unaligned reads analysis\n")
    alignment_rate = open(prefix + "_reads_alignment_rate", 'w')

    num_unaligned = len(unaligned_length)
    if num_unaligned != 0:
        alignment_rate.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
        unaligned_length_2d = unaligned_length[:, numpy.newaxis]
        kde_unaligned = KernelDensity(bandwidth=10).fit(unaligned_length_2d)
        joblib.dump(kde_unaligned, prefix + "_unaligned_length.pkl")
    else:
        alignment_rate.write("Aligned / Unaligned ratio:\t100%\n")

    alignment_rate.close()
    del unaligned_length

    # MATCH AND ERROR MODELS
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
    if args.mode == "transcriptome":
        error_model.hist(prefix + "_genome", alnm_ext)  # Use primary genome alignment for error profiling
    else:
        error_model.hist(prefix, alnm_ext)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        model_fitting.model_fitting(prefix, int(num_threads))

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main()
