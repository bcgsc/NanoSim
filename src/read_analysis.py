#!/usr/bin/env python

"""
@author: Chen Yang, Saber Hafezqorani, Ka Ming Nip, and Theodora Lo
This script generates read profiles for Oxford Nanopore reads (genomic, transcriptomic, and metagenomic).
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
import joblib
from sklearn.neighbors import KernelDensity
import head_align_tail_dist as align
import get_besthit_maf
import get_primary_sam
import besthit_to_histogram as error_model
import model_fitting
import model_intron_retention as model_ir
import model_homopolymer_lengths as model_hp_len
import model_base_qualities as model_base_qual
import pairwise2maf
from file_handler import gzopen as open
from subprocess import DEVNULL


PYTHON_VERSION = sys.version_info
VERSION = "3.2.1"
PROGRAM = "NanoSim"
AUTHOR = "Chen Yang, Saber Hafezqorani, Ka Ming Nip, and Theodora Lo (UBC & BCGSC)"
CONTACT = "cheny@bcgsc.ca; shafezqorani@bcgsc.ca; kmnip@bcgsc.ca"


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
        name, seqs, last = last[1:], [], None
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


def align_transcriptome(in_fasta, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g, chimeric, homopolymer,
                        fastq, quant=False, q_mode=False, normalize=True):
    if t_alnm == '':
        if aligner == "minimap2":
            t_alnm = prefix + "_transcriptome_alnm.bam"
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
            sys.stdout.flush()
            call("minimap2 --cs -ax map-ont -t " + num_threads + " " + ref_t + " " + in_fasta + " | samtools view -b > " + t_alnm,
                 shell=True)

        elif aligner == "LAST":
            t_alnm = prefix + "_transcriptome_alnm.maf"
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference transcriptome\n")
            sys.stdout.flush()
            call("lastdb ref_transcriptome " + ref_t, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_transcriptome " + in_fasta + " > " + t_alnm, shell=True)

    if quant:
        quant = {'tpm': 0}

    # post-process
    t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
    t_alnm_ext = t_alnm_ext[1:]
    sys.stdout.write(
        strftime("%Y-%m-%d %H:%M:%S") + ": Processing transcriptome alignment file: " + t_alnm_ext + '\n')
    sys.stdout.flush()
    if t_alnm_ext == "maf":
        processed_maf_t = prefix + "_transcriptome_alnm_processed.maf"
        call("grep '^s ' " + t_alnm + " > " + processed_maf_t, shell=True)
        unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf_t,
                                                                             prefix + "_transcriptome")
    elif t_alnm_ext in ["sam", "bam"]:
        if chimeric:
            unaligned_length, strandness, unaligned_base_quals = get_primary_sam.primary_and_unaligned_chimeric(t_alnm,
                                                                                        prefix + "_transcriptome",
                                                                                        quant, q_mode, normalize, fastq)
        else:
            unaligned_length, strandness, unaligned_base_quals = get_primary_sam.primary_and_unaligned(t_alnm,
                                                                                                       prefix + "_transcriptome",
                                                                                                       fastq=fastq)

    if q_mode:
        return

    if homopolymer and t_alnm_ext != "maf":
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Converting processed transcriptome alignment file to maf "
                                                         "format for homopolymer length modelling\n")
        sys.stdout.flush()
        call("samtools view " + prefix + "_transcriptome_primary.bam | sam2pairwise > " + prefix + "_transcriptome_primary.out", shell=True, stderr=DEVNULL)
        pairwise2maf.main(["-i" + prefix + "_transcriptome_primary.out", "-o" + prefix + "_transcriptome_alnm_processed.maf"])

    if g_alnm == '':
        if aligner == "minimap2":
            g_alnm = prefix + "_genome_alnm.bam"
            # Alignment to reference genome
            # [EDIT] I may change the options for minimap2 when dealing with cDNA and dRNA reads.
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference genome\n")
            sys.stdout.flush()
            call("minimap2 --cs -ax splice -t " + num_threads + " " + ref_g + " " + in_fasta + " | samtools view -b > " + g_alnm,
                 shell=True)

        elif aligner == "LAST":
            g_alnm = prefix + "_genome_alnm.maf"
            # Alignment to reference genome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference genome\n")
            sys.stdout.flush()
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " > " + g_alnm, shell=True)

    g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
    g_alnm_ext = g_alnm_ext[1:]
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing genome alignment file: " + g_alnm_ext + '\n')
    sys.stdout.flush()
    if g_alnm_ext == "maf":
        processed_maf = prefix + "_processed.maf"
        call("grep '^s ' " + g_alnm + " > " + processed_maf, shell=True)
        get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf, prefix + "_genome")

    elif g_alnm_ext in ["sam", "bam"]:
        get_primary_sam.primary_and_unaligned(g_alnm, prefix + "_genome")

    return t_alnm_ext, unaligned_length, g_alnm, t_alnm, strandness, unaligned_base_quals


def align_genome(in_fasta, prefix, aligner, num_threads, g_alnm, ref_g, chimeric, homopolymer, fastq, quantification=None,
                 q_mode=False):
    # if an alignment file is not provided
    if g_alnm == '':
        if aligner == "minimap2":  # Align with minimap2 by default
            g_alnm = prefix + "_genome_alnm.bam"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2\n")
            sys.stdout.flush()
            call("minimap2 --cs -ax map-ont -t " + num_threads + " " + ref_g + " " + in_fasta + " | samtools view -b > " + g_alnm,
                 shell=True)

        elif aligner == "LAST":
            g_alnm = prefix + "_genome_alnm.maf"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
            sys.stdout.flush()
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " " + g_alnm, shell=True)

    # post-process
    pre, file_ext = os.path.splitext(g_alnm)
    file_extension = file_ext[1:]
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing alignment file: " + file_extension + "\n")
    sys.stdout.flush()
    if file_extension == "maf":
        processed_maf = prefix + "_processed.maf"
        call("grep '^s ' " + g_alnm + " > " + processed_maf, shell=True)

        # get best hit and unaligned reads
        unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf, prefix)

    elif file_extension in ["sam", "bam"]:
        # get the primary alignments and define unaligned reads.
        if chimeric:
            unaligned_length, strandness, unaligned_base_quals = get_primary_sam.primary_and_unaligned_chimeric(g_alnm, prefix,
                                                                                          quantification, q_mode=q_mode, fastq=fastq)
        else:
            unaligned_length, strandness, unaligned_base_quals = get_primary_sam.primary_and_unaligned(g_alnm, prefix, quantification, fastq=fastq)

    if homopolymer and file_extension != "maf":
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Converting processed alignment file to maf format for homopolymer length modelling\n")
        sys.stdout.flush()
        call("samtools view " + prefix + "_primary.bam | sam2pairwise > " + prefix + "_primary.out", shell=True, stderr=DEVNULL)
        pairwise2maf.main(["-i" + prefix + "_primary.out", "-o" + prefix + "_processed.maf"])
    return file_extension, unaligned_length, strandness, unaligned_base_quals


def add_intron(annot, prefix):
    # Read the annotation GTF/GFF3 file
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
    sys.stdout.flush()
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


def concatenate_genomes(g_list):
    metagenome = "reference_metagenome.fasta"
    with open(metagenome, 'wt') as f:
        for g, info in g_list.items():
            path = info['path']
            with open(path, 'rt') as g_f:
                chr = 0
                for line in g_f:
                    if line[0] == '>':
                        f.write('>' + g + '_chr' + str(chr) +'\n')
                        chr += 1
                    else:
                        f.write(line)
    return metagenome


def main():
    parser = argparse.ArgumentParser(
        description=dedent('''
        Read characterization step
        -----------------------------------------------------------
        Given raw ONT reads, reference genome, metagenome, and/or 
        transcriptome, learn read features and output error profiles
        '''),
        formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument('-v', '--version', action='version', version='NanoSim ' + VERSION)
    subparsers = parser.add_subparsers(dest='mode', description=dedent('''
        There are five modes in read_analysis.
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
    parser_g.add_argument('-ga', '--g_alnm', help='Genome alignment file in SAM/BAM format (optional)', default='')
    parser_g.add_argument('-o', '--output', help='The location and prefix of outputting profiles (Default = training)',
                          default='training')
    parser_g.add_argument('-c', '--chimeric', help='Detect chimeric and split reads (Default = False)',
                          action='store_true', default=False)
    parser_g.add_argument('-hp', '--homopolymer', help='Analyze homopolymer lengths (Default = False)',
                          action='store_true', default=False)
    parser_g.add_argument('--min_homopolymer_len', help='Minimum length of homopolymers to analyze (Default = 5 bp)',
                          default='5')
    parser_g.add_argument('--fastq', help='Analyze base qualities (Default = False)', action='store_true', default=False)
    parser_g.add_argument('--no_model_fit', help='Disable model fitting step', action='store_false', default=True)
    parser_g.add_argument('-t', '--num_threads', help='Number of threads for alignment and model fitting (Default = 1)',
                          default='1')

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode")
    parser_t.add_argument('-i', '--read', help='Input read for training', required=True)
    parser_t.add_argument('-rg', '--ref_g', help='Reference genome', required=True)
    parser_t.add_argument('-rt', '--ref_t', help='Reference Transcriptome', required=True)
    parser_t.add_argument('-annot', '--annotation', help='Annotation file in ensemble GTF/GFF formats, '
                                                         'required for intron retention detection', default='')
    parser_t.add_argument('-a', '--aligner', help='The aligner to be used: minimap2 or LAST (Default = minimap2)',
                          choices=['minimap2', 'LAST'], default='minimap2')
    parser_t.add_argument('-ga', '--g_alnm', help='Genome alignment file in SAM/BAM format (optional)', default='')
    parser_t.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in SAM/BAM format (optional)',
                          default='')
    parser_t.add_argument('-o', '--output', help='The location and prefix of outputting profiles (Default = training)',
                          default='training')
    parser_t.add_argument('--no_model_fit', help='Disable model fitting step', action='store_false', default=True)
    parser_t.add_argument('--no_intron_retention', help='Disable Intron Retention analysis', action='store_false',
                          default=True)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads for alignment and model fitting (Default = 1)',
                          default='1')
    parser_t.add_argument('-c', '--chimeric', help='Detect chimeric and split reads (Default = False)',
                          action='store_true', default=False)
    parser_t.add_argument('-q', '--quantification', help='Perform abundance quantification (Default = False)',
                          action='store_true', default=False)
    parser_t.add_argument('-n', '--normalize', help='Normalize by transcript length', action='store_true',
                          default=False)
    parser_t.add_argument('-hp', '--homopolymer', help='Analyze homopolymer lengths (Default = False',
                          action='store_true', default=False)
    parser_t.add_argument('--min_homopolymer_len', help='Minimum length of homopolymers to analyze (Default = 5 bp)',
                          default='5')
    parser_t.add_argument('--fastq', help='Analyze base qualities (Default = False)', action='store_true',
                          default=False)

    parser_m = subparsers.add_parser('metagenome', help="Run the simulator on metagenome mode")
    parser_m.add_argument('-i', '--read', help='Input read for training', required=True)
    parser_m.add_argument('-gl', '--genome_list', help='Reference metagenome list, tsv file, the first column is '
                                                       'species/strain name, the second column is the reference genome '
                                                       'fasta/fastq file directory, the third column is optional, if '
                                                       'provided, it contains the expected abundance (sum up to 100)',
                          required=True)
    parser_m.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam format, the header of each species '
                                                  'should match the metagenome list provided above (optional)',
                          default='')
    parser_m.add_argument('-o', '--output', help='The location and prefix of outputting profiles (Default = training)',
                          default='training')
    parser_m.add_argument('-c', '--chimeric', help='Detect chimeric and split reads (Default = False)',
                          action='store_true', default=False)
    parser_m.add_argument('-q', '--quantification', help='Perform abundance quantification and compute the deviation '
                                                         'between calculated abundance and expected values '
                                                         '(Default = False)',
                          action='store_true', default=False)
    parser_m.add_argument('-hp', '--homopolymer', help='Analyze homopolymer lengths (Default = False)',
                          action='store_true', default=False)
    parser_m.add_argument('--min_homopolymer_len', help='Minimum length of homopolymers to analyze (Default = 5 bp)',
                          default='5')
    parser_m.add_argument('--fastq', help='Analyze base qualities (Default = False)', action='store_true',
                          default=False)
    parser_m.add_argument('--no_model_fit', help='Disable model fitting step', action='store_false', default=True)
    parser_m.add_argument('-t', '--num_threads', help='Number of threads for alignment and model fitting (Default = 1)',
                          default='1')

    parser_e = subparsers.add_parser('quantify', help="Quantify transcriptome expression or metagenome abundance")
    parser_e.add_argument('-e', help='Select from (trans, meta)', default=None, required=True)
    parser_e.add_argument('-i', '--read', help='Input reads for quantification', required=True)
    parser_e.add_argument('-rt', '--ref_t', help='Reference Transcriptome', default='')
    parser_e.add_argument('-gl', '--genome_list', help='Reference metagenome list, tsv file, the first column is '
                                                       'species/strain name, the second column is the reference genome '
                                                       'fasta/fastq file directory, the third column is optional, if '
                                                       'provided, it contains the expected abundance (sum up to 100)',
                          default='')
    parser_e.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam format, the header of each species '
                                                  'should match the metagenome list provided above (optional)',
                          default='')
    parser_e.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in sam format(optional)', default='')
    parser_e.add_argument('-o', '--output', help='The location and prefix of outputting profile (Default = expression)',
                          default='expression')
    parser_e.add_argument('-t', '--num_threads', help='Number of threads for alignment (Default = 1)', default='1')
    parser_e.add_argument('-n', '--normalize', help='Normalize by transcript length', action='store_true',
                          default=False)

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
    parser_ir.add_argument('-ga', '--g_alnm', help='Genome alignment file in SAM/BAM format (optional)', default='')
    parser_ir.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in SAM/BAM format (optional)',
                           default='')
    parser_ir.add_argument('-t', '--num_threads', help='Number of threads for alignment (Default = 1)', default='1')

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    # parse quantify mode arguments
    if args.mode == "quantify":
        infile = args.read
        mode = args.e
        prefix = args.output
        num_threads = args.num_threads
        ref_t = args.ref_t
        genome_list = args.genome_list
        g_alnm = args.g_alnm
        t_alnm = args.t_alnm
        normalize = args.normalize

        print("\nrunning the code with following parameters\n")
        print("mode", mode)
        print("infile", infile)
        print("ref_t", ref_t)
        print("genome_list", genome_list)
        print("g_alnm", g_alnm)
        print("t_aln,", t_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("normalize by transcript length", normalize)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        if mode == "trans":
            # check validity of parameters
            if t_alnm != '':
                pre, file_ext = os.path.splitext(t_alnm)
                file_extension = file_ext[1:]
                if file_extension not in ['sam', 'bam']:
                    print("Please specify an acceptable alignment format! (.sam/.bam)\n")
                    parser_g.print_help(sys.stderr)
                    sys.exit(1)

            # Quantifying the transcript abundance from input read
            align_transcriptome(infile, prefix, 'minimap2', num_threads, t_alnm, ref_t, '', '', True, False, False,
                                True, True, normalize)

        elif mode == "meta":
            # check validity of parameters
            if g_alnm != '':
                pre, file_ext = os.path.splitext(g_alnm)
                file_extension = file_ext[1:]
                if file_extension not in ['sam', 'bam']:
                    print("Please specify an acceptable alignment format! (.sam/.bam)\n")
                    parser_g.print_help(sys.stderr)
                    sys.exit(1)

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing reference genome\n")
            sys.stdout.flush()
            metagenome_list = {}
            with open(genome_list, 'rt') as f:
                for line in f:
                    info = line.strip().split('\t')
                    species = '-'.join(info[0].split())
                    metagenome_list[species] = {'path': info[1]}
                    if len(info) == 3:  # expected abundance
                        metagenome_list[species]['expected'] = float(info[2])
            ref_g = concatenate_genomes(metagenome_list)

            align_genome(infile, prefix, 'minimap2', num_threads, g_alnm, ref_g, True, metagenome_list, True, False)

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
        num_threads = args.num_threads

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
            if file_extension not in ['maf', 'sam', 'bam']:
                print("Please specify an acceptable alignment format! (.maf/.sam /.bam)\n")
                parser_ir.print_help(sys.stderr)
                sys.exit(1)
        if t_alnm != '':
            pre, file_ext = os.path.splitext(t_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam', 'bam']:
                print("Please specify an acceptable alignment format! (.maf/.sam/.bam)\n")
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
        alnm_ext, unaligned_length, g_alnm, t_alnm, strandness, _ \
            = align_transcriptome(infile, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g)

        # Add introns to annotation file
        add_intron(annot, prefix)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
        sys.stdout.flush()
        model_ir.intron_retention(prefix, prefix + "_added_intron_final.gff3", g_alnm, t_alnm)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
        return

    if args.mode == "genome":
        infile = args.read
        ref_g = args.ref_g
        aligner = args.aligner
        g_alnm = args.g_alnm
        prefix = args.output
        num_threads = args.num_threads
        model_fit = args.no_model_fit
        chimeric = args.chimeric
        homopolymer = args.homopolymer
        min_hp_len = args.min_homopolymer_len
        fastq = args.fastq

        # check validity of parameters
        if g_alnm != '':
            pre, file_ext = os.path.splitext(g_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam', 'bam']:
                print("Please specify an acceptable alignment format! (.maf/.sam/.bam)\n")
                parser_g.print_help(sys.stderr)
                sys.exit(1)
        if g_alnm == '' and ref_g == '':
            print("Please supply a reference genome or genome alignment file\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)
        if homopolymer and int(min_hp_len) < 0:
            print("\nPlease input proper min homopolymer length value (>= 0) to analyze\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)
        if fastq:
            reads_pre, reads_file_ext = os.path.splitext(infile)
            if reads_file_ext == ".gz":
                reads_pre, reads_file_ext = os.path.splitext(reads_pre)
            if reads_file_ext not in ['.fastq', '.fq']:
                print("\nFASTA file detected. Please supply reads in FASTQ format (.fastq/.fq/.fastq.gz/.fq.gz).\n")
                sys.exit(1)

        print("\nRunning the code with following parameters:\n")
        print("infile", infile)
        print("ref_g", ref_g)
        print("aligner", aligner)
        print("g_alnm", g_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("model_fit", model_fit)
        print("chimeric", chimeric)
        print("homopolymer", homopolymer)
        if homopolymer:
            print("min_homopolymer_len", min_hp_len)
        print("fastq", fastq)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process\n")
        sys.stdout.flush()
        in_fasta = prefix + "_processed"
        if fastq:
            in_fasta += ".fastq"
        else:
            in_fasta += ".fasta"
        if aligner == 'minimap2':
            in_fasta += '.gz'

        processed_fasta = open(in_fasta, 'wt')

        # Replace spaces in sequence headers with dashes to create unique header for each read
        with open(infile, 'rt') as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                if fastq:
                    processed_fasta.write('@' + chr_name + '\n' + seqS + '\n' + '+' + "\n" + seqQ + '\n')
                else:
                    processed_fasta.write('>' + chr_name + '\n' + seqS + '\n')
        processed_fasta.close()

        alnm_ext, unaligned_length, strandness, unaligned_base_qualities = align_genome(in_fasta, prefix, aligner,
                                                                                        num_threads, g_alnm, ref_g,
                                                                                        chimeric, homopolymer, fastq)
        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        sys.stdout.flush()
        num_aligned = align.head_align_tail(prefix, alnm_ext, args.mode)

    if args.mode == "metagenome":
        infile = args.read
        genome_list = args.genome_list
        g_alnm = args.g_alnm
        quantification = args.quantification
        prefix = args.output
        num_threads = args.num_threads
        model_fit = args.no_model_fit
        chimeric = args.chimeric
        homopolymer = args.homopolymer
        min_hp_len = args.min_homopolymer_len
        fastq = args.fastq

        # check validity of parameters
        if g_alnm != '':
            pre, file_ext = os.path.splitext(g_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['sam', 'bam']:
                print("Please specify an acceptable alignment format! (.sam/.bam)\n")
                parser_g.print_help(sys.stderr)
                sys.exit(1)
        if homopolymer and int(min_hp_len) < 0:
            print("\nPlease input proper min homopolymer length value (>= 0) to analyze\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)
        if fastq:
            reads_pre, reads_file_ext = os.path.splitext(infile)
            if reads_file_ext == ".gz":
                reads_pre, reads_file_ext = os.path.splitext(reads_pre)
            if reads_file_ext not in ['.fastq', '.fq']:
                print("\nFASTA file detected. Please supply reads in FASTQ format (.fastq/.fq/.fastq.gz/.fq.gz).\n")
                sys.exit(1)

        print("\nRunning the code with following parameters:\n")
        print("infile", infile)
        print("genome_list", genome_list)
        print("g_alnm", g_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("model_fit", model_fit)
        print("chimeric", chimeric)
        print("homopolymer", homopolymer)
        if homopolymer:
            print("min_homopolymer_len", min_hp_len)
        print("fastq", fastq)
        print("quantification", quantification)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process\n")
        sys.stdout.flush()
        in_fasta = prefix + "_processed"
        if fastq:
            in_fasta += ".fastq"
        else:
            in_fasta += ".fasta"
        if aligner == 'minimap2':
            in_fasta += '.gz'

        processed_fasta = open(in_fasta, 'wt')
        # Replace spaces in sequence headers with dashes to create unique header for each read
        with open(infile, 'rt') as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                if fastq:
                    processed_fasta.write('@' + chr_name + '\n' + seqS + '\n' + '+' + "\n" + seqQ + '\n')
                else:
                    processed_fasta.write('>' + chr_name + '\n' + seqS + '\n')
        processed_fasta.close()

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing reference genome\n")
        sys.stdout.flush()
        metagenome_list = {}
        with open(genome_list, 'rt') as f:
            for line in f:
                info = line.strip().split('\t')
                species = '-'.join(info[0].split())
                metagenome_list[species] = {'path': info[1]}
                if len(info) == 3:  # expected abundance
                    metagenome_list[species]['expected'] = float(info[2])
        ref_g = concatenate_genomes(metagenome_list)

        metagenome_list = None if not quantification else metagenome_list
        alnm_ext, unaligned_length, strandness, unaligned_base_qualities = align_genome(in_fasta, prefix, 'minimap2',
                                                                                      num_threads, g_alnm, ref_g,
                                                                                      chimeric, homopolymer, fastq,
                                                                                      metagenome_list)
        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        sys.stdout.flush()
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
        num_threads = args.num_threads
        model_fit = args.no_model_fit
        ir = args.no_intron_retention
        chimeric = args.chimeric
        homopolymer = args.homopolymer
        min_hp_len = args.min_homopolymer_len
        fastq = args.fastq
        quantification = args.quantification
        normalize = args.normalize

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
            if g_alnm_ext != t_alnm_ext or g_alnm_ext not in ['maf', 'sam', 'bam']:
                print("\nPlease provide genome and transcriptome alignments in the same format: sam OR bam OR maf\n")
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

        if homopolymer and int(min_hp_len) < 0:
            print("\nPlease input proper min homopolymer length value (>= 0) to analyze\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if fastq:
            reads_pre, reads_file_ext = os.path.splitext(infile)
            print(reads_file_ext)
            if reads_file_ext == ".gz":
                reads_pre, reads_file_ext = os.path.splitext(reads_pre)
            if reads_file_ext not in ['.fastq', '.fq']:
                print("\nFASTA file detected. Please supply reads in FASTQ format (.fastq/.fq/.fastq.gz/.fq.gz).\n")
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
        print("chimeric", chimeric)
        print("homopolymer", homopolymer)
        if homopolymer:
            print("min_homopolymer_len", min_hp_len)
        print("fastq", fastq)
        print("quantification", quantification)
        print("normalize by transcript length", normalize)

        dir_name = os.path.dirname(prefix)
        if dir_name != '':
            call("mkdir -p " + dir_name, shell=True)

        # READ PRE-PROCESS AND ALIGNMENT ANALYSIS
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")
        sys.stdout.flush()
        in_fasta = prefix + "_processed"
        if fastq:
            in_fasta += ".fastq"
        else:
            in_fasta += ".fasta"
        if aligner == 'minimap2':
            in_fasta += '.gz'

        processed_fasta = open(in_fasta, 'wt')
        with open(infile, 'rt') as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                if fastq:
                    processed_fasta.write('@' + chr_name + '\n' + seqS + '\n' + '+' + "\n" + seqQ + '\n')
                else:
                    processed_fasta.write('>' + chr_name + '\n' + seqS + '\n')
        processed_fasta.close()

        alnm_ext, unaligned_length, g_alnm, t_alnm, strandness, unaligned_base_qualities = \
            align_transcriptome(in_fasta, prefix, aligner, num_threads, t_alnm, ref_t, g_alnm, ref_g, chimeric,
                                homopolymer, fastq, min_hp_len, quantification, normalize)

        if ir:
            # Add introns to annotation file
            add_intron(annot, prefix)

            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
            sys.stdout.flush()
            model_ir.intron_retention(prefix, prefix + "_added_intron_final.gff3", g_alnm, t_alnm)

        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        sys.stdout.flush()
        num_aligned = align.head_align_tail(prefix + "_transcriptome", alnm_ext, args.mode)

    # strandness of the aligned reads
    strandness_rate = open(prefix + "_strandness_rate", 'wt')
    strandness_rate.write("strandness:\t" + str(round(strandness, 3)))
    strandness_rate.close()

    # Length distribution of unaligned reads
    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Unaligned reads analysis\n")
    sys.stdout.flush()
    alignment_rate = open(prefix + "_reads_alignment_rate", 'wt')

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
    sys.stdout.flush()
    if args.mode == "transcriptome":
        error_model.hist(prefix + "_genome", alnm_ext)  # Use primary genome alignment for error profiling
    else:
        error_model.hist(prefix, alnm_ext)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        sys.stdout.flush()
        model_fitting.model_fitting(prefix, int(num_threads))

    # HOMOPOLYMER LENGTH MODEL
    if homopolymer:
        sys.stdout.write(
            strftime("%Y-%m-%d %H:%M:%S") + ": Analyzing homopolymer lengths and estimating model parameters\n")
        sys.stdout.flush()
        if args.mode == "transcriptome":
            model_hp_len.model_homopolymer_lengths(prefix + "_transcriptome_alnm_processed.maf", min_hp_len, prefix)
        else:
            model_hp_len.model_homopolymer_lengths(prefix + "_processed.maf", min_hp_len, prefix)

    # BASE QUALITIES MODEL
    if fastq:
        if alnm_ext == "maf":  # TODO: implement base quality model parameter estimation from alignments in MAF format
            sys.stderr.write("Currently, base quality model parameter estimation cannot be performed on alignments in "
                             "MAF format.\n")
            sys.exit(1)
        else:
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Analyzing base qualities and estimating model parameters\n")
            sys.stdout.flush()
            if args.mode == "transcriptome":
                model_base_qual.model_base_qualities(prefix + "_transcriptome_primary." + alnm_ext, prefix,
                                                     unaligned_base_qualities)
            else:
                model_base_qual.model_base_qualities(prefix + "_primary." + alnm_ext, prefix, unaligned_base_qualities)

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")

if __name__ == "__main__":
    main()
