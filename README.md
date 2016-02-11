# NanoSim  

## Dependencies
LAST (Tested with version 581)
R (Tested with version 3.2.3)

## Characterization stage  
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

## Simulation stage  
__Usage:__  
```
./simulator.py [command] <options>  
   [command]:  
    circular | linear  
    # Do not choose 'circular' when there is more than one sequence in the reference  
    <options>:  
    -h : print usage message
    -r : reference genome in fasta file, specify path and file name  
    -c : Flowcell model_prefix, same as the output prefix in read_analysis.py, default = training  
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
