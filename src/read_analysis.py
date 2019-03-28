#!/usr/bin/env python

"""
@author: Chen Yang & Saber HafezQorani

This script generates read profiles for Oxford Nanopore 2D reads (genomic and transcriptome).

"""


from __future__ import print_function
from __future__ import with_statement
from subprocess import call
from time import strftime
try:
    from six.moves import xrange
except ImportError:
    pass
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


def align_transcriptome(in_fasta, prefix, aligner, num_threads, g_alnm, t_alnm, ref_t, ref_g):

    if g_alnm != "" and t_alnm != "":
        out_g = g_alnm
        out_t = t_alnm
        g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
        t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
        g_alnm_ext = g_alnm_ext [1:]
        t_alnm_ext = t_alnm_ext[1:]
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Processing the alignment files: " + t_alnm_ext + "\n")
        if g_alnm_ext == t_alnm_ext == "maf":
            processed_maf_g = prefix + "_genome_alnm_processed.maf"
            processed_maf_t = prefix + "_transcriptome_alnm_processed.maf"
            call("grep '^s ' " + g_alnm + " > " + processed_maf_g, shell=True)
            call("grep '^s ' " + t_alnm + " > " + processed_maf_t, shell=True)

            unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, processed_maf_t, prefix)

        elif g_alnm_ext == t_alnm_ext == "sam":

            unaligned_length, strandness = get_primary_sam.primary_and_unaligned(t_alnm, prefix)

    elif (g_alnm == '' and t_alnm == ''):
        if aligner == "minimap2":
            g_alnm_ext = "sam"
            t_alnm_ext = "sam"
            outsam_g = prefix + "_genome_alnm.sam"
            outsam_t = prefix + "_transcriptome_alnm.sam"
            out_g = outsam_g
            out_t = outsam_t

            # Alignment to reference genome
            # [EDIT] I may change the options for minimap2 when dealing with cDNA and dRNA reads.
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference genome\n")
            call("minimap2 -ax splice " + ref_g + " " + in_fasta + " > " + outsam_g, shell=True)
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2 to reference transcriptome\n")
            call("minimap2 --cs -ax map-ont " + ref_t + " " + in_fasta + " > " + outsam_t, shell=True)

            # [EDIT] I may add a script to remove minimap2/LAST post-alignment files after alignment.
            unaligned_length, strandness = get_primary_sam.primary_and_unaligned(outsam_t, prefix)

        elif aligner == "LAST":
            g_alnm_ext = "maf"
            t_alnm_ext = "maf"
            outmaf_g = prefix + "_genome_alnm.maf"
            outmaf_t = prefix + "_transcriptome_alnm.maf"
            out_g = outmaf_g
            out_t = outmaf_t
            # Alignment to reference genome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference genome\n")
            call("lastdb ref_genome " + ref_g, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " | grep '^s ' > " + outmaf_g,
                 shell=True)
            # Alignment to reference transcriptome
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST to reference transcriptome\n")
            call("lastdb ref_transcriptome " + ref_t, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_transcriptome " + in_fasta + " | grep '^s ' > " + outmaf_t,
                 shell=True)

            unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, outmaf_t, prefix)

    return t_alnm_ext, unaligned_length, out_g, out_t, strandness


def align_genome(in_fasta, prefix, aligner, num_threads, g_alnm, ref_g):
    # if an alignment file is provided
    if g_alnm != '':
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


    # if alignment file is not provided
    else:
        if aligner == "minimap2" or aligner == "":  # Align with minimap2 by default
            file_extension = "sam"
            out_sam = prefix + "_genome_alnm.sam"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with minimap2\n")
            call("minimap2 --cs -ax map-ont -t " + num_threads + " " + ref_g + " " + in_fasta + " > " + out_sam,
                 shell=True)
            # get primary alignments and unaligned reads
            unaligned_length, strandness = get_primary_sam.primary_and_unaligned(out_sam, prefix)
        elif aligner == "LAST":
            file_extension = "maf"
            out_maf = prefix + "_genome_alnm.maf"
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
            call("lastdb ref_genome " + ref, shell=True)
            call("lastal -a 1 -P " + num_threads + " ref_genome " + in_fasta + " | grep '^s ' > " + out_maf, shell=True)
            unaligned_length, strandness = get_besthit_maf.besthit_and_unaligned(in_fasta, out_maf, prefix)

    return file_extension, unaligned_length, strandness


def main(argv):
    # Parse input and output files

    prefix = 'training'
    model_fit = True
    intron_retention = True
    detect_IR = False
    quantify = False

    parser = argparse.ArgumentParser(
        description='Given the read profiles from characterization step, ' \
                    'simulate genomic/transcriptic ONT reads and output error profiles',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    subparsers = parser.add_subparsers(help = "You may run the simulator on transcriptome or genome mode. You may also only quanity expression profiles.", dest='mode')

    parser_g = subparsers.add_parser('genome', help="Run the simulator on genome mode.")
    parser_g.add_argument('-i', '--read', help='Input read for training.', required=False)
    parser_g.add_argument('-rg', '--ref_g', help='Reference genome.', required=False)
    parser_g.add_argument('-a', '--aligner', help='The aligner to be used minimap2 or LAST (Default = minimap2)', default = 'minimap2')
    parser_g.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format (optional)', default= '')
    parser_g.add_argument('-o', '--output', help='The output name and location for profiles', default = "training")
    parser_g.add_argument('--no_model_fit', help='Disable model fitting step', action='store_true')
    parser_g.add_argument('-t', '--num_threads', help='Number of threads to be used in alignments and model fitting (Default = 1)', default=1)

    parser_t = subparsers.add_parser('transcriptome', help="Run the simulator on transcriptome mode.")
    parser_t.add_argument('-i', '--read', help='Input read for training.', required=False)
    parser_t.add_argument('-rg', '--ref_g', help='Reference genome.', required=False)
    parser_t.add_argument('-rt', '--ref_t', help='Reference Transcriptome.', required=False)
    parser_t.add_argument('-annot', '--annot', help='Annotation file in ensemble GTF/GFF formats.', required=False)
    parser_t.add_argument('-a', '--aligner', help='The aligner to be used minimap2 or LAST (Default = minimap2)', default = 'minimap2')
    parser_t.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format (optional)', default= '')
    parser_t.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in sam or maf format (optional)', default= '')
    parser_t.add_argument('-o', '--output', help='The output name and location for profiles', default = "training")
    parser_t.add_argument('--no_model_fit', help='Disable model fitting step', action='store_true')
    parser_t.add_argument('--no_intron_retention', help='Disable Intron Retention analysis', action='store_true')
    parser_t.add_argument('-b', '--num_bins', help='Number of bins to be used (Default = 20)', default = 20)
    parser_t.add_argument('-t', '--num_threads', help='Number of threads to be used in alignments and model fitting (Default = 1)', default=1)

    parser_e = subparsers.add_parser('quantify', help="Quantify expression profile of transcripts")
    parser_e.add_argument('-o', '--output', help='The output name and location', default="training")
    parser_e.add_argument('-i', '--read', help='Input reads to use for quantification.', required=True)
    parser_e.add_argument('-rt', '--ref_t', help='Reference Transcriptome.', required=True)
    parser_e.add_argument('-t', '--num_threads', help='Number of threads to be used (Default = 1)', default=1)

    parser_ir = subparsers.add_parser('detect_ir', help="Detect Intron Retention events using the alignment file")
    parser_ir.add_argument('-annot', '--annot', help='Annotation file in ensemble GTF/GFF formats.', required=True)
    parser_ir.add_argument('-o', '--output', help='The output name and location', default = "ir_info", required=True)
    parser_ir.add_argument('-ga', '--g_alnm', help='Genome alignment file in sam or maf format', default= '', required=True)
    parser_ir.add_argument('-ta', '--t_alnm', help='Transcriptome alignment file in sam or maf format', default= '', required=True)

    args = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    if len(sys.argv) == 2:
        if args.mode == "genome":
            parser_g.print_help(sys.stderr)
        elif args.mode == "transcriptome":
            parser_t.print_help(sys.stderr)
        elif args.mode == "detect_ir":
            parser_ir.print_help(sys.stderr)
        elif args.mode == "quantify":
            parser_e.print_help(sys.stderr)
        else:
            parser.print_help(sys.stderr)
        sys.exit(1)


    #parse quanity mode arguments
    if args.mode == "quantify":
        infile = args.read
        ref_t = args.ref_t
        prefix = args.output
        num_threads = max(args.num_threads, 1)

        print("running the code with following parameters:\n")
        print("infile", infile)
        print("ref_t", ref_t)
        print("prefix", prefix)
        print("num_threads", num_threads)

        # Quantifying the transcript abundance from input read
        sys.stdout.write('Quantifying transcripts abundance: \n')
        call("minimap2 -t " + str(num_threads) + " -x map-ont -p0 " + ref_t + " " + infile + " > " + prefix + "_mapping.paf", shell=True)
        call("python nanopore_transcript_abundance.py -i " + prefix + "_mapping.paf > " + prefix + "_abundance.tsv",
             shell=True)
        sys.stdout.write('Finished! \n')
        sys.exit(1)

    if args.mode == "detect_ir":
        annot = args.annot
        prefix = args.output
        g_alnm = args.g_alnm
        t_alnm = args.t_alnm

        if g_alnm == "" or t_alnm == "":
            print("Please provide both alignments in sam format\n")
            parser_ir.print_help(sys.stderr)
            sys.exit(1)

        print("running the code with following parameters:\n")
        print("annot", annot)
        print("g_alnm", g_alnm)
        print("t_alnm", t_alnm)
        print("prefix", prefix)

        # Read the annotation GTF/GFF3 file
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
        # If gtf provided, convert to GFF3 (gt gtf_to_gff3)
        annot_filename, annot_file_extension = os.path.splitext(annot)
        annot_file_extension = annot_file_extension[1:]
        if annot_file_extension.upper() == "GTF":
            call("gt gtf_to_gff3 -tidy -o " + prefix + ".gff3 " + annot, shell=True)

        # Next, add intron info into gff3:
        call("gt gff3 -tidy -retainids -checkids -addintrons -force -o " + prefix + "_addedintron.gff3 " + annot_filename + ".gff3",
            shell=True)
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
        model_ir.intron_retention(prefix, prefix + "_addedintron.gff3", g_alnm, t_alnm)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")
        sys.exit(1)

    if args.mode == "genome":
        infile = args.read
        ref_g = args.ref_g
        aligner = args.aligner
        g_alnm = args.g_alnm
        prefix = args.output
        num_threads = max(args.num_threads, 1)
        if args.no_model_fit:
            model_fit = False

        if aligner not in ['minimap2', 'LAST', '']:
            print("Please specify an acceptable aligner (minimap2 or LAST)\n")
            parser_g.print_help(sys.stderr)
            sys.exit(1)

        if g_alnm != '':
            pre, file_ext = os.path.splitext(g_alnm)
            file_extension = file_ext[1:]
            if file_extension not in ['maf', 'sam']:
                print("Please specify an acceptable alignment format! (.maf or .sam)\n")
                parser_g.print_help(sys.stderr)
                sys.exit(1)

        print("running the code with following parameters:\n")
        print("infile", infile)
        print("ref_g", ref_g)
        print("aligner", aligner)
        print("g_alnm", g_alnm)
        print("prefix", prefix)
        print("num_threads", num_threads)
        print("model_fit", model_fit)

        dir_name = os.path.dirname(prefix)
        basename = os.path.basename(prefix)
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

        alnm_ext, unaligned_length, strandness = align_genome(in_fasta, prefix, aligner, num_threads, g_alnm, ref_g)

        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        num_aligned = align.head_align_tail(prefix, alnm_ext, args.mode)

    if args.mode == "transcriptome":
        infile = args.read
        ref_g = args.ref_g
        ref_t = args.ref_t
        annot = args.annot
        aligner = args.aligner
        g_alnm = args.g_alnm
        t_alnm = args.t_alnm
        prefix = args.output
        num_bins = max(args.num_bins, 20) #I may remove it because of ecdf > KDE
        num_threads = max(args.num_threads, 1)
        if args.no_model_fit:
            model_fit = False
        if args.no_intron_retention:
            intron_retention = False

        if aligner not in ['minimap2', 'LAST', '']:
            print("Please specify an acceptable aligner (minimap2 or LAST)\n")
            parser_t.print_help(sys.stderr)
            sys.exit(1)

        if (g_alnm != '' and t_alnm == '') or (g_alnm == '' and t_alnm != ''):
            print("Please specify either both alignment files (-ga and -ta) OR an aligner to use for alignment (-a)")
            parser_t.print_help(sys.stderr)
            sys.exit(1)
        if g_alnm != "" and t_alnm != "":
            g_alnm_filename, g_alnm_ext = os.path.splitext(g_alnm)
            t_alnm_filename, t_alnm_ext = os.path.splitext(t_alnm)
            g_alnm_ext = g_alnm_ext[1:]
            t_alnm_ext = t_alnm_ext[1:]
            if g_alnm_ext != t_alnm_ext:
                print("Please provide both alignments in a same format: sam OR maf\n")
                parser_t.print_help(sys.stderr)
                sys.exit(1)


        print("running the code with following parameters:\n")
        print("infile", infile)
        print("ref_g", ref_g)
        print("ref_t", ref_t)
        print("annot", annot)
        print("aligner", aligner)
        print("g_alnm", g_alnm)
        print("t_alnm", t_alnm)
        print("prefix", prefix)
        print("num_bins", num_bins)
        print("num_threads", num_threads)
        print("model_fit", model_fit)
        print("intron_retention", intron_retention)

        dir_name = os.path.dirname(prefix)
        basename = os.path.basename(prefix)
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

        alnm_ext, unaligned_length, out_g, out_t, strandness = align_transcriptome(in_fasta, prefix, aligner, num_threads, g_alnm, t_alnm, ref_t, ref_g)

        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read the length of reference transcripts \n")
        # Read the length of reference transcripts from the reference transcriptome
        dict_ref_len = {}
        with open(ref_t) as f:
            for seqN, seqS, seqQ in readfq(f):
                info = re.split(r'[_\s]\s*', seqN)
                chr_name = "-".join(info)
                dict_ref_len[chr_name] = len(seqS)

        if intron_retention:
            # Read the annotation GTF/GFF3 file
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Parse the annotation file (GTF/GFF3)\n")
            # If gtf provided, convert to GFF3 (gt gtf_to_gff3)
            annot_filename, annot_file_extension = os.path.splitext(annot)
            annot_file_extension = annot_file_extension[1:]
            if annot_file_extension.upper() == "GTF":
                call("gt gtf_to_gff3 -tidy -o " + prefix + ".gff3 " + annot, shell=True)

            # Next, add intron info into gff3:
            call("gt gff3 -tidy -retainids -checkids -addintrons -force -o " + prefix + "_addedintron.gff3 " + annot_filename + ".gff3",
                shell=True)
            sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Modeling Intron Retention\n")
            model_ir.intron_retention(prefix, prefix + "_addedintron.gff3", out_g, out_t)

        # Aligned reads analysis
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
        num_aligned = align.head_align_tail(prefix, alnm_ext, args.mode, dict_ref_len)

    # strandness of the aligned reads
    strandness_rate = open(prefix + "_strandness_rate", 'w')
    strandness_rate.write("strandness:\t" + str(round(strandness, 3)))
    strandness_rate.close()

    # Length distribution of unaligned reads
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
    error_model.hist(prefix, alnm_ext)

    if model_fit:
        sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
        model_fitting.model_fitting(prefix, int(num_threads))

    sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
    main(sys.argv[1:])

