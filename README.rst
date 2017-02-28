NanoSimH
========

.. image:: https://travis-ci.org/karel-brinda/NanoSimH.svg?branch=master
		:target: https://travis-ci.org/karel-brinda/NanoSimH

NanoSimH is a modified version of `NanoSim`_, a fast and scalable read simulator that captures the technology-specific features of ONT data, and allows for adjustments upon improvement of nanopore sequencing technology. It has been created as a fork of NanoSim 1.0.1. The main improvements compared to `NanoSim`_ are the following:

* Support for Python 3
* Support for `RNF`_ read names
* Installation from `PyPI`_
* Automatic testing using `Travis`_
* Reproducible simulations (setting a seed for PRG)
* Improved interface with new parameters (e.g., for merging all contigs)
* Several minor bugs fixed

.. _RNF: https://www.ncbi.nlm.nih.gov/pubmed/26353839
.. _PyPI: https://pypi.python.org/pypi/NanoSimH/
.. _Travis: https://travis-ci.org/karel-brinda/NanoSimH
.. _NanoSim: https://github.com/bcgsc/NanoSim

Installation
------------

**From PyPI**

.. code-block:: bash

	pip install --upgrade nanosimh

**From git**

.. code-block:: bash

		git clone https://github.com/karel-brinda/nanosimh
		cd nanosimh
		pip install --upgrade .

or

.. code-block:: bash

		git clone https://github.com/karel-brinda/nanosimh
		cd nanosimh
		python setup.py install


Dependencies
------------

* LAST (Tested with version 581)  
* R (Tested with version 3.2.3)  
* Python (2.7, 3.2 - 3.6)  
* Numpy (Tested with version 1.10.1 or above)  

Usage
-----

NanoSim is implemented using R for error model fitting and Python for read length analysis and simulation. The first step of NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference genome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

1. Characterization stage
~~~~~~~~~~~~~~~~~~~~~~~~~

Characterization stage takes a reference and a training read set in FASTA format as input. User can also provide their own alignment file in MAF format.  

**Usage:**


.. code-block::

	$ nanosimh_train --help
	usage: nanosimh_train [-h] [-i str] -r str [-m str] [-o str] [-b int]
	                      [--no-model-fit]

	NanoSimH - a fork of NanoSim, a simulator of Oxford Nanopore reads.

	optional arguments:
	  -h, --help            show this help message and exit
	  -i str, --infile str  training ONT real reads, must be fasta files
	  -r str, --ref str     reference genome of the training reads
	  -m str, --maf str     user can provide their own alignment file, with maf
	                        extension
	  -o str                prefix of output file [training]
	  -b int                number of bins (for development) [20]
	  --no-model-fit        no model fitting


\* NOTICE: -m option allows users to provide their own alignment file. Make sure that the name of query sequences are the same as appears in the fasta files. For fasta files, some headers have spaces in them and most aligners only take part of the header (before the first white space/tab) as the query name. However, the truncated headers may not be unique if using the output of poretools. We suggest users to pre-process the fasta files by concatenating all elements in the header via '\_' before alignment and feed the processed fasta file as input of NanoSim.  

Some ONT read profiles are ready to use for users. With the profiles, users can run simulation tool directly. Please go to ftp://ftp.bcgsc.ca/supplementary/NanoSim/ to download *E. coli* or *S. cerevisiae* datasets and profiles.

2. Simulation stage  
~~~~~~~~~~~~~~~~~~~

Simulation stage takes reference genome and read profiles as input and outputs simulated reads in FASTA fomat.  

**Usage:**

.. code-block::

	$ nanosimh_simulate --help
	usage: nanosimh_simulate [-h] -r str [-p str] [-o str] [-n int] [-m float]
	                         [-i float] [-d float] [-s int] [--circular]
	                         [--perfect] [--merge-contigs] [--rnf]
	                         [--rnf-add-cigar] [--max-len int] [--min-len int]
	                         [--kmer-bias int]

	NanoSimH - a fork of NanoSim, a simulator of Oxford Nanopore reads.

	optional arguments:
	  -h, --help            show this help message and exit
	  -r str, --reference str
	                        reference genome in fasta file
	  -p str, --profile str
	                        prefix of training set profiles [training]
	  -o str, --out-pref str
	                        prefix of output file [simulated]
	  -n int, --number int  number of generated reads [20000]
	  -m float, --mis-rate float
	                        mismatch rate (weight tuning) [1.0]
	  -i float, --ins-rate float
	                        insertion rate (weight tuning) [1.0]
	  -d float, --del-rate float
	                        deletion reate (weight tuning) [1.0]
	  -s int, --seed int    initial seed for the pseudorandom number generator (0
	                        for random) [1]
	  --circular            circular simulation (linear otherwise)
	  --perfect             output perfect reads, no mutations
	  --merge-contigs       merge contigs from the reference
	  --rnf                 use RNF format for read names
	  --rnf-add-cigar       add cigar to RNF names
	  --max-len int         maximum read length [inf]
	  --min-len int         minimum read length [50]
	  --kmer-bias int       prohibits homopolymers with length >= n bases in
	                        output reads [6]

	Notice: the use of `max_len` and `min_len` will affect the read length
	distributions. If the range between `max_len` and `min_len` is too small, the
	program will run slowlier accordingly.  

**For example:**

1 If you want to simulate *E. coli* genome, then circular command must be chosen because it's a circular genome  
``nanosimh_simulate --circular -r Ecoli_ref.fasta -p ecoli``

2 If you want to simulate only perfect reads, i.e. no snps, or indels, just simulate the read length distribution  
``nanosimh_simulate --circular -r Ecoli_ref.fasta -p ecoli --perfect``

3 If you want to simulate *S. cerevisiae* genome with no kmer bias, then linear command must be chosen because it's a linear genome  
``nanosimh_simulate -r yeast_ref.fasta -p yeast --kmer-bias 0``

*See more detailed example in example.sh*

Explaination of output files  
----------------------------

1. Characterization stage
~~~~~~~~~~~~~~~~~~~~~~~~~

1. ``training_aligned_length_ecdf`` Length distribution of aligned regions on aligned reads  
2. ``training_aligned_reads_ecdf`` Length distribution of aligned reads  
3. ``training_align_ratio`` Empirical distribution of align ratio of each read  
4. ``training_besthit.maf`` The best alignment of each read based on length  
5. ``training_match.hist/training_mis.hist/training_del.hist/training_ins.hist`` Histogram of match, mismatch, and indels  
6. ``training_first_match.hist`` Histogram of the first match length of each alignment  
7. ``training_error_markov_model`` Markov model of error types  
8. ``training_ht_ratio`` Empirical distribution of the head region vs total unaligned region  
9. ``training.maf`` The output of LAST, alignment file in MAF format  
10. ``training_match_markov_model`` Markov model of the length of matches (stretches of correct base calls)  
11. ``training_model_profile`` Fitted model for errors  
12. ``training_processed.maf`` A re-formatted MAF file for user-provided alignment file  
13. ``training_unaligned_length_ecdf`` Length distribution of unaligned reads  

2. Simulation stage  
~~~~~~~~~~~~~~~~~~~

1. ``simulated.log``

	Log file for simulation process  
	
2. ``simulated_reads.fasta``

	FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
	
	To explain the information in the header, we have two examples:  

	* ``>ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0``  
		All information before the first ``_`` are chromosome information. ``468529`` is the start position and *unaligned* suggesting it should be unaligned to the reference. The first ``0`` is the sequence index. ``F`` represents a forward strand. ``0_3236_0`` means that sequence length extracted from the reference is 3236 bases.  
	* ``>ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2``
		This is an aligned read coming from chromosome XI at position 115406. ``16565`` is the sequence index. `R` represents a reverse complement strand. ``92_12710_2`` means that this read has 92-base head region (cannot be aligned), followed by 12710 bases of middle region, and then 2-base tail region.  
	
	The information in the header can help users to locate the read easily.  
	
3. ``simulated_error_profile``

	Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases.  
