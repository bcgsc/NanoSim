1. Characterization stage
Usage:
./read_analysis.py <options>  
    [options]:  
    -h : print usage message  
    -i : training ONT real reads, must be fasta files  
    -r : reference genome of the training reads  
    -o : The prefix of output file, default = 'training'

2. Simulation stage
Usage:
./simulator.py [command] <options> 
    [command]:
    circular | linear
    # Do not choose 'circular' when there is more than one sequence in the reference
    <options>: 
    -r: reference genome in fasta file, specify path and file name
    -c : Flowcell chemistry, R7 or R7.3
    -o : The prefix of output file, default = 'simulated'
    -n : Number of generated reads, default = 24,221 reads
    -s : Substitution profile, can be omitted if there is no customized profile
    -p : Error model profile, can be omitted if there is no customized profile
    -i : Insertion rate, a floating number in the interval [0, 1], default = 0.05
    -d : Deletion rate, a floating number in the interval [0, 1], default = 0.05
    -m : Mismatch rate, a floating number in the interval [0, 1], default = 0.1
    --perfect: Output perfect reads, no mutations, default = False

For example:
1. If you want to simulate E. coli genome, then circular command must be chosen because it's a circular genome
./simulator.py circular -r Ecoli_ref.fasta -c R7 

2. If you want to simulate only perfect reads, i.e.no snps, or indels, just simulate the read length distribution
./simulator.py circular -r Ecoli_ref.fasta -c R7 --perfect
