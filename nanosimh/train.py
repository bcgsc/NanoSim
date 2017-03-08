#!/usr/bin/env python

"""
Created on Apr 10, 2015

@author: Chen Yang

This script generates read profiles Oxford Nanopore 2D reads.

"""

from __future__ import print_function
from __future__ import with_statement
from subprocess import call
from time import strftime
try:
	test_xrange=xrange(42)
except NameError:
	from six.moves import xrange
import sys
import os
import argparse
import numpy
from .head_align_tail_dist import *
from .get_besthit import *
from .besthit_to_histogram import * 
import multiprocessing
from .misc import *

nb_cores=multiprocessing.cpu_count()

def run(command):
	print("Running '{}'".format(command),file=sys.stderr)
	call(command, shell=True)

def main():
	# Parse input and output files
	infile = ''
	outfile = 'training'
	ref = ''
	maf_file = ''
	model_fit = True
	num_bins = 20

	parser = argparse.ArgumentParser(
			description='NanoSimH - a fork of NanoSim, a simulator of Oxford Nanopore reads.',
		)

	parser.add_argument('-i', '--infile',
			type=str,
			metavar='str',
			dest='infile',
			help='training ONT real reads, must be fasta files',
			default='',
		)
	parser.add_argument('-r', '--ref',
			type=str,
			metavar='str',
			required=True,
			dest='ref',
			help='reference genome of the training reads',
		)
	parser.add_argument('-m', '--maf',
			type=str,
			metavar='str',
			dest='maf_file',
			help='user can provide their own alignment file, with maf extension',
			default='',
		)
	parser.add_argument('-p','--profile',
			type=str,
			metavar='str',
			dest='outfile',
			help='prefix of output files [{}]'.format(outfile),
			default=outfile,
		)
	parser.add_argument('-b','--num-bins',
			type=int,
			metavar='int',
			dest='num_bins',
			help='number of bins (for development) [{}]'.format(num_bins),
			default=num_bins,
		)
	parser.add_argument('--no-model-fit',
			action='store_false',
			dest='model_fit',
			help='no model fitting',
		)

	args = parser.parse_args()

	infile = args.infile
	outfile = args.outfile
	ref = args.ref
	maf_file = args.maf_file
	model_fit = args.model_fit
	num_bins = args.num_bins

	assert num_bins>0

	assert infile!='' or ref!=''
	if infile!='':
		assert_file_exists(infile, True)
	if ref!='':
		assert_file_exists(ref, True)


	# READ PRE-PROCESS AND UNALIGNED READS ANALYSIS
	sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Read pre-process and unaligned reads analysis\n")
	out_maf = outfile + ".maf"

	# Read pre-process
	in_fasta = outfile + ".fasta"
	if in_fasta == infile:
		in_fasta = outfile + "_processed.fasta"
	out_fasta = open(in_fasta, 'w')
	dic_reads = {}
	with open(infile, 'r') as f:
		for line in f:
			if line[0] == '>':
				name = '-'.join(line.strip()[1:].split())
				dic_reads[name] = ""
			else:
				dic_reads[name] += line.strip()
	for k, v in dic_reads.items():
		out_fasta.write('>' + k + '\n' + v + '\n')
	out_fasta.close()

	del dic_reads

	# if maf file provided
	if maf_file != '':
		assert_file_exists(maf_file, True)
		if out_maf == maf_file:
			out_maf = outfile + "_processed.maf"

		call("grep '^s ' \"{}\" > \"{}\"".format(maf_file, out_maf), shell=True)

		# get best hit and unaligned reads
		unaligned_length = list(besthit_and_unaligned(in_fasta, out_maf, outfile))

	# if maf file not provided
	else:
		# Alignment
		sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Alignment with LAST\n")
		run('lastdb -P {} ref_genome "{}"'.format(nb_cores, ref))
		run('lastal -T0 -r1 -q1 -a1 -b1 -m100 -P {} ref_genome "{}" | grep \'^s \' > "{}"'.format(nb_cores, in_fasta, out_maf))

		# get best hit and unaligned reads
		unaligned_length = list(besthit_and_unaligned(in_fasta, out_maf, outfile))

	# ALIGNED READS ANALYSIS
	sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Aligned reads analysis\n")
	num_aligned = head_align_tail(outfile, num_bins)

	# Length distribution of unaligned reads
	out1 = open(outfile + "_unaligned_length_ecdf", 'w')

	num_unaligned = len(unaligned_length)
	if num_unaligned != 0:
		max_length = max(unaligned_length)
		hist_unaligned, edges_unaligned = numpy.histogram(unaligned_length, bins=numpy.arange(0, max_length + 50, 50),
														  density=True)
		cdf = numpy.cumsum(hist_unaligned * 50)
		out1.write("Aligned / Unaligned ratio:" + "\t" + str(num_aligned * 1.0 / num_unaligned) + '\n')
		out1.write("bin\t0-" + str(max_length) + '\n')
		for i in xrange(len(cdf)):
			out1.write(str(edges_unaligned[i]) + '-' + str(edges_unaligned[i+1]) + "\t" + str(cdf[i]) + '\n')
	else:
		out1.write("Aligned / Unaligned ratio:\t100%\n")

	out1.close()
	del unaligned_length

	# MATCH AND ERROR MODELS
	sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": match and error models\n")
	hist(outfile)

	if model_fit:
		sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Model fitting\n")
		r_path = os.path.join(os.path.dirname(__file__), "model_fitting.R")
		if os.path.isfile(r_path):
			run("R CMD BATCH '--args prefix=\"{}\"' \"{}\"".format(outfile,r_path))
		else:
			sys.stderr.write("Could not find 'model_fitting.R' in ../src/\n" +
				  "Make sure you copied the whole source files from Github.")

	sys.stdout.write(strftime("%Y-%m-%d %H:%M:%S") + ": Finished!\n")


if __name__ == "__main__":
	main()

