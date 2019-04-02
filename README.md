![NanoSim](https://github.com/bcgsc/NanoSim/blob/master/NanoSim%20logo.png)

NanoSim is a fast and scalable read simulator that captures the technology-specific features of ONT data, and allows for adjustments upon improvement of nanopore sequencing technology.  

The second version of NanoSim (v2.0) uses minimap2 as default aligner to align long genomic ONT reads to reference genome. It leads to much faster alignment step and reduces the overall runtime of NanoSim. We also utilize HTSeq, a python package, to read SAM alignment files efficiently.

The latest version of NanoSim (v2.5) is able to simulate ONT transcriptome reads (cDNA / directRNA) as well as genomic reads. It also models features of the library preparation protocols used, including intron retention (IR) events in cDNA and directRNA reads. Further, it optionally profiles transcript expression patterns.

__Citation__: Chen Yang, Justin Chu, René L Warren, Inanç Birol; NanoSim: nanopore sequence read simulator based on statistical characterization. Gigascience 2017 gix010. doi: 10.1093/gigascience/gix010

## Dependencies
minimap2 (Tested with version 2.10)  
LAST (Tested with version 581 and 916)  
Python (2.7 or >= 3.4)  
Python packages:  
* six  
* numpy (Tested with version 1.10.1 or above)
* HTSeq  
* scipy (Tested with verson 1.0.0)  
* scikit-learn (Tested with version 0.20.0)

## Usage
NanoSim is implemented using Python for error model fitting, read length analysis, and simulation. The first step of NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference genome/transcriptome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

### 1. Characterization stage  
Characterization stage runs in four mode: genome, transcriptome, quantify, and detect_ir. Below you may see the general usage of this code. We will explain each mode separately as well.

__Characterization step general usage:__
```
usage: read_analysis.py [-h] {genome,transcriptome,quantify,detect_ir} ...

Given the read profiles from characterization step, simulate
genomic/transcriptic ONT reads and output error profiles

positional arguments:
  {genome,transcriptome,quantify,detect_ir}
                        You may run the simulator on transcriptome or genome
                        mode. You may also only quanity expression profiles.
    genome              Run the simulator on genome mode.
    transcriptome       Run the simulator on transcriptome mode.
    quantify            Quantify expression profile of transcripts
    detect_ir           Detect Intron Retention events using the alignment
                        file

optional arguments:
  -h, --help            show this help message and exit
```


If you are interested in simulating ONT genomic reads, you need to run the characterization stage in "genome" mode with following options. It takes a reference genome and a training read set in FASTA format as input and aligns these reads to the reference using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats. The output of this is a bunch of profiles which you should use in simulation stage.
__genome mode usage:__
```
usage: read_analysis.py genome [-h] -i READ -rg REF_G [-a ALIGNER]
                               [-ga G_ALNM] [-o OUTPUT] [--no_model_fit]
                               [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i READ, --read READ  Input read for training.
  -rg REF_G, --ref_g REF_G
                        Reference genome.
  -a ALIGNER, --aligner ALIGNER
                        The aligner to be used minimap2 or LAST (Default =
                        minimap2)
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam or maf format (optional)
  -o OUTPUT, --output OUTPUT
                        The output name and location for profiles
  --no_model_fit        Disable model fitting step
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads to be used in alignments and model
                        fitting (Default = 1)
```

If you are interested in simulating ONT transcriptome reads (cDNA / directRNA), you need to run the characterization stage in "transcriptome" mode with following options. It takes a reference transcriptome and a training read set in FASTA format as input and aligns these reads to the reference using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats. The output of this is a bunch of profiles which you should use in simulation stage.
__transcriptome mode usage:__
```
usage: read_analysis.py transcriptome [-h] -i READ [-rg REF_G] -rt REF_T
                                      [-annot ANNOT] [-a ALIGNER] [-ga G_ALNM]
                                      [-ta T_ALNM] [-o OUTPUT]
                                      [--no_model_fit] [--no_intron_retention]
                                      [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i READ, --read READ  Input read for training.
  -rg REF_G, --ref_g REF_G
                        Reference genome.
  -rt REF_T, --ref_t REF_T
                        Reference Transcriptome.
  -annot ANNOT, --annot ANNOT
                        Annotation file in ensemble GTF/GFF formats.
  -a ALIGNER, --aligner ALIGNER
                        The aligner to be used: minimap2 or LAST (Default =
                        minimap2)
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam or maf format (optional)
  -ta T_ALNM, --t_alnm T_ALNM
                        Transcriptome alignment file in sam or maf format
                        (optional)
  -o OUTPUT, --output OUTPUT
                        The output name and location for profiles
  --no_model_fit        Disable model fitting step
  --no_intron_retention
                        Disable Intron Retention analysis
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads to be used in alignments and model
                        fitting (Default = 1)
```

The "transcriptome" mode of the NanoSim is able to model features of the library preparation protocols used, including intron retention (IR) events in cDNA and directRNA reads. Further, it optionally profiles transcript expression patterns. However, if you are interested in only detecting Intron Retention events or quantifying expression patterns of transcripts without running other characterization stage, you may use two modes we introduced for this purpose: "quantify" and "detect_ir". Details are as follows:

__quantifty mode usage:__
```
usage: read_analysis.py quantify [-h] [-o OUTPUT] -i READ -rt REF_T
                                 [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -o OUTPUT, --output OUTPUT
                        The output name and location
  -i READ, --read READ  Input reads to use for quantification.
  -rt REF_T, --ref_t REF_T
                        Reference Transcriptome.
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads to be used (Default = 1)
```

__detect_ir mode usage:__
```
usage: read_analysis.py detect_ir [-h] -annot ANNOT [-o OUTPUT] -ga G_ALNM -ta
                                  T_ALNM

optional arguments:
  -h, --help            show this help message and exit
  -annot ANNOT, --annot ANNOT
                        Annotation file in ensemble GTF/GFF formats.
  -o OUTPUT, --output OUTPUT
                        The output name and location
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam or maf format
  -ta T_ALNM, --t_alnm T_ALNM
                        Transcriptome alignment file in sam or maf format
```


\* NOTICE: -ga/-ta option allows users to provide their own alignment file. Make sure that the name of query sequences are the same as appears in the FASTA files. For FASTA files, some headers have spaces in them and most aligners only take part of the header (before the first white space/tab) as the query name. However, the truncated headers may not be unique if using the output of poretools. We suggest users to pre-process the fasta files by concatenating all elements in the header via '\_' before alignment and feed the processed FASTA file as input of NanoSim.

~~Some ONT read profiles are ready to use for users. With the profiles, users can run simulation tool directly. Please go to [ftp](ftp://ftp.bcgsc.ca/supplementary/NanoSim/) to download _E. coli_ or _S. cerevisiae_ datasets and profiles.~~  
UPDATE: For release v2.2.0 onwards, there's no pre-trained profiles to download temporarily. We are working on this for now.

### 2. Simulation stage
Simulation stage takes reference genome/transcriptome and read profiles as input and outputs simulated reads in FASTA format. Simulation stage runs in two modes: "genome" and "transcriptome" and you may use either of them based on your needs.

__Simulationb stage general usage:__
```
usage: simulator.py [-h] {genome,transcriptome} ...

Given the read profiles from characterization step, simulate
genomeic/transcriptomic ONT reads and outputs the error profiles.

positional arguments:
  {genome,transcriptome}
                        You may run the simulator on transcriptome or genome
                        mode.
    genome              Run the simulator on genome mode.
    transcriptome       Run the simulator on transcriptome mode.

optional arguments:
  -h, --help            show this help message and exit
```
If you are interested in simulating ONT genomic reads, you need to run the simulation stage in "genome" mode with following options.
__genome mode usage:__
```
usage: simulator.py genome [-h] -rg REF_G [-c MODEL_PREFIX] [-o OUTPUT]
                           [-n NUMBER] [-i INSERTION_RATE] [-d DELETION_RATE]
                           [-m MISMATCH_RATE] [-max MAX_LEN] [-min MIN_LEN]
                           [-med MEDIAN_LEN] [-sd SD_LEN] [--seed SEED]
                           [-k KMERBIAS] [-s STRANDNESS] [-dna_type DNA_TYPE]
                           [--perfect]

optional arguments:
  -h, --help            show this help message and exit
  -rg REF_G, --ref_g REF_G
                        Input reference genome
  -c MODEL_PREFIX, --model_prefix MODEL_PREFIX
                        Address for profiles created in characterization step
                        (model_prefix)
  -o OUTPUT, --output OUTPUT
                        Output address for simulated reads
  -n NUMBER, --number NUMBER
                        Number of reads to be simulated
  -i INSERTION_RATE, --insertion_rate INSERTION_RATE
                        Insertion rate (optional)
  -d DELETION_RATE, --deletion_rate DELETION_RATE
                        Deletion rate (optional)
  -m MISMATCH_RATE, --mismatch_rate MISMATCH_RATE
                        Mismatch rate (optional)
  -max MAX_LEN, --max_len MAX_LEN
                        The maximum length for simulated reads
  -min MIN_LEN, --min_len MIN_LEN
                        The minimum length for simulated reads
  -med MEDIAN_LEN, --median_len MEDIAN_LEN
                        The median read length
  -sd SD_LEN, --sd_len SD_LEN
                        The standard deviation of read length in log scale
  --seed SEED           Manually seeds the pseudo-random number generator
  -k KMERBIAS, --KmerBias KMERBIAS
                        Determine whether to considert Kmer Bias or not
  -s STRANDNESS, --strandness STRANDNESS
                        Determine the strandness of the simulated reads.
                        Overrides the value profiled in characterization
                        phase. Should be between 0 and 1.
  -dna_type DNA_TYPE    Specify the dna type: circular OR linear, default =
                        linear
  --perfect             Ignore profiles and simulate perfect reads
```


If you are interested in simulating ONT transcriptome reads, you need to run the simulation stage in "transcriptome" mode with following options.
__transcriptome mode usage:__
```
usage: simulator.py transcriptome [-h] -rt REF_T -rg REF_G [-e EXP]
                                  [-c MODEL_PREFIX] [-o OUTPUT] [-n NUMBER]
                                  [-i INSERTION_RATE] [-d DELETION_RATE]
                                  [-m MISMATCH_RATE] [-max MAX_LEN]
                                  [-min MIN_LEN] [-k KMERBIAS] [-s STRANDNESS]
                                  [--no_model_ir] [--perfect]

optional arguments:
  -h, --help            show this help message and exit
  -rt REF_T, --ref_t REF_T
                        Input reference transcriptome
  -rg REF_G, --ref_g REF_G
                        Input reference genome
  -e EXP, --exp EXP     Expression profile in the specified format specified
                        in the documentation
  -c MODEL_PREFIX, --model_prefix MODEL_PREFIX
                        Address for profiles created in characterization step
                        (model_prefix)
  -o OUTPUT, --output OUTPUT
                        Output address for simulated reads
  -n NUMBER, --number NUMBER
                        Number of reads to be simulated
  -i INSERTION_RATE, --insertion_rate INSERTION_RATE
                        Insertion rate (optional)
  -d DELETION_RATE, --deletion_rate DELETION_RATE
                        Deletion rate (optional)
  -m MISMATCH_RATE, --mismatch_rate MISMATCH_RATE
                        Mismatch rate (optional)
  -max MAX_LEN, --max_len MAX_LEN
                        The maximum length for simulated reads
  -min MIN_LEN, --min_len MIN_LEN
                        The minimum length for simulated reads
  -k KMERBIAS, --KmerBias KMERBIAS
                        Determine whether to considert Kmer Bias or not
  -s STRANDNESS, --strandness STRANDNESS
                        Determine the strandness of the simulated reads.
                        Overrides the value profiled in characterization
                        phase. Should be between 0 and 1.
  --no_model_ir         Consider Intron Retention model from characterization
                        step when simulating reads
  --perfect             Ignore profiles and simulate perfect reads
```


\* Notice: the use of `max_len` and `min_len` will affect the read length distributions. If the range between `max_len` and `min_len` is too small, the program will run slowlier accordingly.

__For example:__  
1 If you want to simulate _E. coli_ genome, then circular command must be chosen because it's a circular genome  
`./simulator.py genome -dna_type circular -rg Ecoli_ref.fasta -c ecoli`

2 If you want to simulate only perfect reads, _i.e._ no snps, or indels, just simulate the read length distribution  
`./simulator.py genome -dna_type circular -rg Ecoli_ref.fasta -c ecoli --perfect`

3 If you want to simulate _S. cerevisiae_ genome with kmer bias, then linear command must be chosen because it's a linear genome  
`./simulator.py genome -dna_type linear -rg yeast_ref.fasta -c yeast --KmerBias`

_See more detailed example in example.sh_  

## Explanation of output files
### 1. Characterization stage
#### 1.1 Characterization stage (genome)
1. `training_aligned_region.pkl` Kernel density function of aligned regions on aligned reads  
2. `training_aligned_reads.pkl` Kernel density function of aligned reads  
3. `training_ht_length.pkl`  Kernel density function of unaligned regions on aligned reads  
4. `training_besthit.maf/sam` The best alignment of each read based on length  
5. `training_match.hist/training_mis.hist/training_del.hist/training_ins.hist` Histogram of match, mismatch, and indels  
6. `training_first_match.hist` Histogram of the first match length of each alignment  
7. `training_error_markov_model` Markov model of error types  
8. `training_ht_ratio.pkl` Kernel density function of head/(head + tail) on aligned reads    
9. `training.maf/sam` The alignment output  
10. `training_match_markov_model` Markov model of the length of matches (stretches of correct base calls)  
11. `training_model_profile` Fitted model for errors  
12. `training_processed.maf` A re-formatted MAF file for user-provided alignment file  
13. `training_unaligned_length.pkl` Kernel density function of unaligned reads  
14. `training_error_rate.tsv` Mismatch rate, insertion rate and deletion rate
15. `training_strandness_rate` Strandness rate in input reads.

#### 1.1 Characterization stage (transcriptome)
1. `training_aligned_region.pkl` Kernel density function of aligned regions on aligned reads
2. `training_aligned_region_2d.pkl` Two-dimensional kernel density function of aligned regions over the length of reference transcript they aligned
3. `training_aligned_reads.pkl` Kernel density function of aligned reads
4. `training_ht_length.pkl`  Kernel density function of unaligned regions on aligned reads
5. `training_besthit.maf/sam` The best alignment of each read based on length
6. `training_match.hist/training_mis.hist/training_del.hist/training_ins.hist` Histogram of match, mismatch, and indels
7. `training_first_match.hist` Histogram of the first match length of each alignment
8. `training_error_markov_model` Markov model of error types
9. `training_ht_ratio.pkl` Kernel density function of head/(head + tail) on aligned reads
10. `training.maf/sam` The alignment output
11. `training_match_markov_model` Markov model of the length of matches (stretches of correct base calls)
12. `training_model_profile` Fitted model for errors
13. `training_processed.maf` A re-formatted MAF file for user-provided alignment file
14. `training_unaligned_length.pkl` Kernel density function of unaligned reads
15. `training_error_rate.tsv` Mismatch rate, insertion rate and deletion rate
16. `training_strandness_rate` Strandness rate in input reads.
17. `training_addedintron_final.gff3` gff3 file format containing the intron coordinate information
18. `training_IR_info` List of transcripts in which there is a retained intron based on IR modeling step
19. `training_IR_markov_model` Markov model of Intron Retention events


### 2. Simulation stage  

1. `simulated_reads.fasta`
  FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
  
  To explain the information in the header, we have two examples:  
  * `>ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0`  
    All information before the first `_` are chromosome information. `468529` is the start position and `unaligned` suggesting it should be unaligned to the reference. The first `0` is the sequence index. `F` represents a forward strand. `0_3236_0` means that sequence length extracted from the reference is 3236 bases.  
  * `>ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2`  
    This is an aligned read coming from chromosome XI at position 115406. `16565` is the sequence index. `R` represents a reverse complement strand. `92_12710_2` means that this read has 92-base head region (cannot be aligned), followed by 12710 bases of middle region, and then 2-base tail region.  
  
  The information in the header can help users to locate the read easily.  
  
2. `simulated_error_profile`
  Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases.  


## Acknowledgements
Sincere thanks to our labmates and all contributors and users of this tool.
