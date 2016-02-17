# NanoSim  
NanoSim is a fast and scalable read simulator that captures the technology-specific features of ONT data, and allows for adjustments upon improvement of nanopore sequencing technology.

## Dependencies
LAST (Tested with version 581)  
R (Tested with version 3.2.3)  
Python (2.6 or above)  
Numpy (Tested with version 1.10.1 or above)  

## Usage
NanoSim is implemented using R for error model fitting and Python for read length analysis and simulation. The first step of NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference genome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

### 1. Characterization stage  
Characterization stage takes a reference and a training read set in FASTA format as input. User can also provide their own alignment file in MAF format.  

__Usage:__  
```
./read_analysis.py <options>  
    [options]:  
    -h : print usage message  
    -i : training ONT real reads, must be fasta files  
    -r : reference genome of the training reads  
    -m : User can provide their own alignment file, with maf extension, can be omitted  
    -o : The prefix of output file, default = 'training'  
```
Some ONT read profiles are ready to use for users. With the profiles, users can run simulation tool directly. Please go to [ftp://ftp.bcgsc.ca/supplementary/NanoSim/] to download _E. coli_ or _S. cerevisiae_ datasets and profiles.

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
    --perfect: Output perfect reads, no mutations, default = False  
    --KmerBias: prohibits homopolymers with length >= 6 bases in output reads, can be omitted  
```
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
1. `training_aligned_length_ecdf` Length distribution of aligned regions on aligned reads  
2. `training_aligned_reads_ecdf` Length distribution of aligned reads  
3. `training_align_ratio` Empirical distribution of align ratio of each read  
4. `training_besthit.maf` The best alignment of each read based on length  
5. `training_match.hist/training_mis.hist/training_del.hist/training_ins.hist` Histogram of match, mismatch, and indels  
6. `training_first_match.hist` Histogram of the first match length of each alignment  
7. `training_error_markov_model` Markov model of error types  
8. `training_ht_ratio` Empirical distribution of the head region vs total unaligned region  
9. `training.maf` The output of LAST, alignment file in MAF format  
10. `training_match_markov_model` Markov model of the length of matches (stretches of correct base calls)  
11. `training_model_profile` Fitted model for errors  
12. `training_processed.maf` A re-formatted MAF file for user-provided alignment file  
13. `training_unaligned_length_ecdf` Length distribution of unaligned reads  

### 2. Simulation stage  
1. `simulated.log`  
  Log file for simulation process  
  
2. `simulated_reads.fasta`  
  FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
  
  To explain the information in the header, we have two examples:  
  * `>ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0`  
    All information before the first `_` are chromosome information. `468529` is the start position and `unaligned` suggesting it should be unaligned to the reference. The first `0` is the sequence index. `F` represents a forward strand. `0_3236_0` means that sequence length extracted from the reference is 3236 bases.  
  * `>ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2`  
    This is an aligned read coming from chromosome XI at position 115406. `16565` is the index of simulation. `R` represents a reverse complement strand. `92_12710_2` means that this read has 92-base head region (cannot be aligned), followed by 12710 bases of middle region, and then 2-base tail region.  
  
  The information in the header can help users to locate the read easily.  
  
3. `simulated_error_profile`  
  Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases.  
