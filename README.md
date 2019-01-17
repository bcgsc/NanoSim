![NanoSim](https://github.com/bcgsc/NanoSim/blob/master/NanoSim%20logo.png)

NanoSim is a fast and scalable read simulator that captures the technology-specific features of ONT data, and allows for adjustments upon improvement of nanopore sequencing technology.  

The second version of NanoSim (v2.0) uses minimap2 as default aligner to align long genomic ONT reads to reference genome. It leads to much faster alignment step and reduces the overall runtime of NanoSim. We also utilize HTSeq, a python package, to read SAM alignment files efficiently.

__Citation__: Chen Yang, Justin Chu, René L Warren, Inanç Birol; NanoSim: nanopore sequence read simulator based on statistical characterization. Gigascience 2017 gix010. doi: 10.1093/gigascience/gix010

## Transcriptome simulation
For transcriptome simulaton (directRNA / cDNA reads), please use [Trans-NanoSim](https://github.com/bcgsc/Trans-NanoSim) for now. We are working on merging Trans-NanoSim pipeline into NanoSim repository and it will be available shortly.

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
NanoSim is implemented using Python for error model fitting, read length analysis, and simulation. The first step of NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference genome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

### 1. Characterization stage  
Characterization stage takes a reference and a training read set in FASTA format as input and aligns these reads to the reference using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats.  

__Usage:__  
```
./read_analysis.py <options>  
    [options]:  
    -h : print usage message  
    -i : training ONT real reads, must be fasta files  
    -r : reference genome of the training reads  
    -a : Aligner to be used: minimap2 or LAST, default = 'minimap2' 
    -m : User can provide their own alignment file, with maf or sam extension, can be omitted  
    -o : The prefix of output file, default = 'training'  
    -t : number of threads for alignment and model fitting, default = 1  
    --no_model_fit : Skip the model fitting step  
```
\* NOTICE: -m option allows users to provide their own alignment file. Make sure that the name of query sequences are the same as appears in the fasta files. For fasta files, some headers have spaces in them and most aligners only take part of the header (before the first white space/tab) as the query name. However, the truncated headers may not be unique if using the output of poretools. We suggest users to pre-process the fasta files by concatenating all elements in the header via '\_' before alignment and feed the processed fasta file as input of NanoSim.  

~~Some ONT read profiles are ready to use for users. With the profiles, users can run simulation tool directly. Please go to [ftp](ftp://ftp.bcgsc.ca/supplementary/NanoSim/) to download _E. coli_ or _S. cerevisiae_ datasets and profiles.~~  
UPDATE: For release v2.2.0 onwards, there's no pre-trained profiles to download temporarily. We are working on this for now.

### 2. Simulation stage  
Simulation stage takes reference genome and read profiles as input and outputs simulated reads in FASTA fomat.  

__Usage:__  
```
./simulator.py [command] <options>  
   [command]:  
    circular | linear  
    # Do not choose 'circular' when there is more than one sequence in the reference  
    <options>:  
    -h : print usage message
    -r : reference genome in fasta file, specify path and file name  
    -c : the prefix of training set profiles, same as the output prefix in read_analysis.py, default = training  
    -o : The prefix of output file, default = 'simulated'  
    -n : Number of generated reads, default = 20,000 reads  
    --max_len : Maximum read length, default = Inf
    --min_len : Minimum read length, default = 50
    --perfect: Output perfect reads, no mutations, default = False  
    --KmerBias: prohibits homopolymers with length >= 6 bases in output reads, can be omitted  
    --seed : manually seeds the pseudo-random number generator, default = None  
    --median_len : the median read length, default = None  
    --sd_len : the standard deviation of read length in log scale, default = None  
```
\* Notice: the use of `max_len` and `min_len` will affect the read length distributions. If the range between `max_len` and `min_len` is too small, the program will run slowlier accordingly.  

__For example:__  
1 If you want to simulate _E. coli_ genome, then circular command must be chosen because it's a circular genome  
`./simulator.py circular -r Ecoli_ref.fasta -c ecoli`  

2 If you want to simulate only perfect reads, _i.e._ no snps, or indels, just simulate the read length distribution  
`./simulator.py circular -r Ecoli_ref.fasta -c ecoli --perfect` 

3 If you want to simulate _S. cerevisiae_ genome with kmer bias, then linear command must be chosen because it's a linear genome  
`./simulator.py linear -r yeast_ref.fasta -c yeast --KmerBias`  

_See more detailed example in example.sh_  

## Explaination of output files  
### 1. Characterization stage
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

### 2. Simulation stage  
1. `simulated.log`  
  Log file for simulation process  
  
2. `simulated_reads.fasta`  
  FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
  
  To explain the information in the header, we have two examples:  
  * `>ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0`  
    All information before the first `_` are chromosome information. `468529` is the start position and `unaligned` suggesting it should be unaligned to the reference. The first `0` is the sequence index. `F` represents a forward strand. `0_3236_0` means that sequence length extracted from the reference is 3236 bases.  
  * `>ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2`  
    This is an aligned read coming from chromosome XI at position 115406. `16565` is the sequence index. `R` represents a reverse complement strand. `92_12710_2` means that this read has 92-base head region (cannot be aligned), followed by 12710 bases of middle region, and then 2-base tail region.  
  
  The information in the header can help users to locate the read easily.  
  
3. `simulated_error_profile`  
  Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases.  


## Acknowledgements
Sincere thanks to my labmates and all contributors and users to this tool.
