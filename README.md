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
