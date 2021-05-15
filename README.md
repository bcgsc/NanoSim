[![Release](https://img.shields.io/github/v/release/bcgsc/nanosim?include_prereleases)](https://github.com/bcgsc/NanoSim/releases)
[![Downloads](https://img.shields.io/github/downloads/bcgsc/Nanosim/total?logo=github)](https://github.com/bcgsc/NanoSim/archive/v3.0.0.zip)
[![Conda](https://img.shields.io/conda/dn/bioconda/nanosim?label=Conda)](https://anaconda.org/bioconda/nanosim)
[![Stars](https://img.shields.io/github/stars/bcgsc/NanoSim.svg)](https://github.com/bcgsc/NanoSim/stargazers)  

![NanoSim](https://github.com/bcgsc/NanoSim/blob/master/NanoSim_logo.png)

NanoSim is a fast and scalable read simulator that captures the technology-specific features of ONT data, and allows for adjustments upon improvement of nanopore sequencing technology.  

The second version of NanoSim (v2.0.0) uses minimap2 as default aligner to align long genomic ONT reads to reference genome. It leads to much faster alignment step and reduces the overall runtime of NanoSim. We also utilize HTSeq, a python package, to read SAM alignment files efficiently.

NanoSim [(v2.5)](https://github.com/bcgsc/NanoSim/releases/tag/v2.5.1) is able to simulate ONT transcriptome reads (cDNA / direct RNA) as well as genomic reads. It also models features of the library preparation protocols used, including intron retention (IR) events in cDNA and directRNA reads. Further, it has stand-alone modes which profiles transcript expression patterns and detects IR events in custom datasets. Additionally, we improved the homopolymer simulation option which simulates homopolymer expansion and contraction events with respect to chosen basecaller. Multiprocessing option allows for faster runtime for large library simulation.  

NanoSim [(v2.6)](https://github.com/bcgsc/NanoSim/releases/tag/v2.6.0) is able to simulate ONT reads in fastq format. The base quality information is simulated with truncated log-normal distributions, learnt separately from match bases, mismatch bases, insertion bases, deletion bases, and unaligned bases, each from different basecaller and read type.  

NanoSim [(v3.0)](https://github.com/bcgsc/NanoSim/releases/tag/v3.0.0) is able to simulate ONT metagenome reads. It quantifies metagenome abundance in the characterization stage, and accomodates for chimeric reads. In the simulation stage, it simulates both features as well. In addition, the simulation of chimeric reads is available in genome mode too. Some pre-trained models are re-trained for compatibility issues.

**We provide 9 pre-trained models in the latest release! Users can choose to download the whole package or only scripts without models to speed it up**

![Citation](https://img.shields.io/badge/NanoSim-manuscript-ff69b4)  
If you use NanoSim to simulate genomic reads, please cite the following paper:

**NanoSim**  
Chen Yang, Justin Chu, René L Warren, and Inanç Birol; NanoSim: nanopore sequence read simulator based on statistical characterization. GigaScience, Volume 6, Issue 4, April 2017, gix010, https://doi.org/10.1093/gigascience/gix010


If you use NanoSim to simulate transcriptomic reads, please cite the following paper:

**Trans-NanoSim**  
Saber Hafezqorani, Chen Yang, Theodora Lo, Ka Ming Nip, René L. Warren, and Inanç Birol; Trans-NanoSim characterizes and simulates nanopore RNA-sequencing data. GigaScience, Volume 9, Issue 6, June 2020, giaa061, https://doi.org/10.1093/gigascience/giaa061


## Dependencies
![Python}](https://img.shields.io/pypi/pyversions/py)  
Python packages:  
* six  
* numpy (Tested with version 1.10.1 or above)
* HTSeq (Tested with version 0.9.1)  
* Pysam (Tested with version 0.13)  
* scipy (Tested with verson 1.0.0)  
* scikit-learn (Tested with version 0.20.0)

minimap2 (Tested with version 2.10 and 2.17)  
LAST (Tested with version 581 and 916)  
samtools (Tested with version 1.12)

## Usage
NanoSim is implemented using Python for error model fitting, read length analysis, and simulation. The first step of NanoSim is read characterization, which provides a comprehensive alignment-based analysis, and generates a set of read profiles serving as the input to the next step, the simulation stage. The simulation tool uses the model built in the previous step to produce in silico reads for a given reference genome/transcriptome. It also outputs a list of introduced errors, consisting of the position on each read, error type and reference bases.

### 1. Characterization stage  
Characterization stage runs in five mode: genome, transcriptome, metagenome, quantify, and detect_ir. Below you may see the general usage of this code. We will explain each mode separately as well.

__Characterization step general usage:__
```
usage: read_analysis.py [-h] [-v]
                        {genome,transcriptome,metagenome, quantify,detect_ir} ...

Read characterization step
-----------------------------------------------------------
Given raw ONT reads, reference genome and/or transcriptome,
learn read features and output error profiles

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  
  There are five modes in read_analysis.
  For detailed usage of each mode:
      read_analysis.py mode -h
  -------------------------------------------------------

  {genome,transcriptome,quantify,detect_ir}
    genome              Run the simulator on genome mode
    transcriptome       Run the simulator on transcriptome mode
    metagenome          Run the simulator on metagenome mode
    quantify            Quantify transcriptome expression or metagenome
                        abundance
    detect_ir           Detect Intron Retention events using the alignment
                        file

```

**genome mode**  
If you are interested in simulating ONT genomic reads, you need to run the characterization stage in "genome" mode with following options. It takes a reference genome and a training read set in FASTA or FASTQ format as input and aligns these reads to the reference using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats. If the SAM file is provided, make sure that is MD flag in the SAM file. The output of this is a bunch of profiles which you should use in simulation stage.

__genome mode usage:__
```
usage: read_analysis.py genome [-h] -i READ [-rg REF_G] [-a {minimap2,LAST}]
                               [-ga G_ALNM] [-o OUTPUT] [-c] [--no_model_fit]
                               [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i READ, --read READ  Input read for training
  -rg REF_G, --ref_g REF_G
                        Reference genome, not required if genome alignment
                        file is provided
  -a {minimap2,LAST}, --aligner {minimap2,LAST}
                        The aligner to be used, minimap2 or LAST (Default =
                        minimap2)
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam/bam or maf format (optional)
  -o OUTPUT, --output OUTPUT
                        The location and prefix of outputting profiles
                        (Default = training)
  -c, --chimeric        Detect chimeric and split reads (Default = False)
  --no_model_fit        Disable model fitting step
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for alignment and model fitting
                        (Default = 1)

```

**transcriptome mode**  
If you are interested in simulating ONT transcriptome reads (cDNA / directRNA), you need to run the characterization stage in "transcriptome" mode with following options. It takes a reference transcriptome, a reference genome, and a training read set in FASTA or FASTQ format as input and aligns these reads to the reference using minimap2 (default) or LAST aligner. User can also provide their own alignment file in SAM or MAF formats. If the SAM file is provided, make sure that is MD flag in the SAM file. The output of this is a bunch of profiles which you should use in simulation stage.

__transcriptome mode usage:__
```
usage: read_analysis.py transcriptome [-h] -i READ [-rg REF_G] -rt REF_T
                                      [-annot ANNOTATION] [-a {minimap2,LAST}]
                                      [-ga G_ALNM] [-ta T_ALNM] [-o OUTPUT]
                                      [--no_model_fit] [--no_intron_retention]
                                      [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i READ, --read READ  Input read for training
  -rg REF_G, --ref_g REF_G
                        Reference genome
  -rt REF_T, --ref_t REF_T
                        Reference Transcriptome
  -annot ANNOTATION, --annotation ANNOTATION
                        Annotation file in ensemble GTF/GFF formats, required
                        for intron retention detection
  -a {minimap2,LAST}, --aligner {minimap2,LAST}
                        The aligner to be used: minimap2 or LAST (Default =
                        minimap2)
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam/bam or maf format (optional)
  -ta T_ALNM, --t_alnm T_ALNM
                        Transcriptome alignment file in sam/bam or maf format
                        (optional)
  -o OUTPUT, --output OUTPUT
                        The location and prefix of outputting profiles
                        (Default = training)
  --no_model_fit        Disable model fitting step
  --no_intron_retention
                        Disable Intron Retention analysis
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for alignment and model fitting
                        (Default = 1)

```

**metagenome mode**
If you are interested in simulating ONT metagenome reads, you need to run the characterization stage in "metagenome" mode with following options. It takes a metagenome list with paths pointing to each genome and a training read set in FASTA or FASTQ format as input and aligns these reads to the reference using minimap2. User can also provide their own alignment file in SAM formats. If the SAM file is provided, make sure that is MD flag in the SAM file. The output of this is a bunch of profiles which you should use in simulation stage.

__metagenome mode usage:__
```
usage: read_analysis.py metagenome [-h] -i READ -gl GENOME_LIST [-ga G_ALNM]
                                   [-o OUTPUT] [-c] [-q] [--no_model_fit]
                                   [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -i READ, --read READ  Input read for training
  -gl GENOME_LIST, --genome_list GENOME_LIST
                        Reference metagenome list, tsv file, the first column
                        is species/strain name, the second column is the
                        reference genome fasta/fastq file directory, the third
                        column is optional, if provided, it contains the
                        expected abundance (sum up to 100)
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam/bam format, the header of
                        each species should match the metagenome list provided
                        above (optional)
  -o OUTPUT, --output OUTPUT
                        The location and prefix of outputting profiles
                        (Default = training)
  -c, --chimeric        Detect chimeric and split reads (Default = False)
  -q, --quantification  Perform Salmon quantification and compute the
                        variation in abundance when compared to expected
                        values (Default = False)
  --no_model_fit        Disable model fitting step
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for alignment and model fitting
                        (Default = 1)
```


**quantify mode**  
If you are interested in quantifying the expression levels of transcriptome or abundance of metagenome, you can run this mode with following options. This mode is independent from the other features in characterization stage, and the output is not used for simulation stage.

__quantifty mode usage:__
```
usage: read_analysis.py quantify [-h] -e E -i READ [-rt REF_T]
                                 [-gl GENOME_LIST] [-ga G_ALNM] [-o OUTPUT]
                                 [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -e E                  Select from (trans, meta)
  -i READ, --read READ  Input reads for quantification
  -rt REF_T, --ref_t REF_T
                        Reference Transcriptome
  -gl GENOME_LIST, --genome_list GENOME_LIST
                        Reference metagenome list, tsv file, the first column
                        is species/strain name, the second column is the
                        reference genome fasta/fastq file directory, the third
                        column is optional, if provided, it contains the
                        expected abundance (sum up to 100)
  -a G_ALNM, --alnm G_ALNM
                        Genome alignment file in sam format, the header of
                        each species should match the metagenome list provided
                        above (optional)
  -o OUTPUT, --output OUTPUT
                        The location and prefix of outputting profile (Default
                        = expression)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for alignment (Default = 1)

```

If `-e trans` is used, users need to provide the reference transcriptome through `-rt` parameter; if `-e meta` is used, users need to provide a genome list containing all reference genome and the path to them through `-gl` parameter OR provide genome alignment file through `-a` parameter.

**detect_ir mode usage**  
The "transcriptome" mode of the NanoSim is able to model features of the library preparation protocols used, including intron retention (IR) events in cDNA and directRNA reads. Further, it optionally profiles transcript expression patterns. However, if you are interested in only detecting Intron Retention events without running other analysis in the characterization stage, you may use the "detect_ir" mode. Details are as follows:

__detect_ir mode usage:__
```
usage: read_analysis.py detect_ir [-h] -annot ANNOTATION [-i READ] [-rg REF_G]
                                  [-rt REF_T] [-a {minimap2,LAST}] [-o OUTPUT]
                                  [-ga G_ALNM] [-ta T_ALNM] [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -annot ANNOTATION, --annotation ANNOTATION
                        Annotation file in ensemble GTF/GFF formats
  -i READ, --read READ  Input read for training, not required if alignment
                        files are provided
  -rg REF_G, --ref_g REF_G
                        Reference genome, not required if genome alignment
                        file is provided
                        
  -rt REF_T, --ref_t REF_T
                        Reference Transcriptome, not required if transcriptome
                        alignment file is provided
  -a {minimap2,LAST}, --aligner {minimap2,LAST}
                        The aligner to be used: minimap2 or LAST (Default =
                        minimap2)
  -o OUTPUT, --output OUTPUT
                        The output name and location
  -ga G_ALNM, --g_alnm G_ALNM
                        Genome alignment file in sam/bam or maf format (optional)
  -ta T_ALNM, --t_alnm T_ALNM
                        Transcriptome alignment file in sam/bam or maf format
                        (optional)
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for alignment (Default = 1)

```

\* NOTICE: -ga/-ta option allows users to provide their own alignment file. Make sure that the name of query sequences are the same as appears in the FASTA files. For FASTA files, some headers have spaces in them and most aligners only take part of the header (before the first white space/tab) as the query name. However, the truncated headers may not be unique if using the output of poretools. We suggest users to pre-process the fasta files by concatenating all elements in the header via '\_' before alignment and feed the processed FASTA file as input of NanoSim.

### Downloads  
**Some ONT read profiles are ready to use for users. With the profiles, users can run simulation tool directly.**  

For **releases before v2.2.0**, we provide profiles trained for _E. coli_ or _S. cerevisiae_ datasets. Flowcell chemistry is R7.3 and R9, and they were basecalled by Metrichor. They can be downloaded from **[our ftp site](http://www.bcgsc.ca/downloads/supplementary)**

For **release v2.5.0 and onwards**, we provide profiles trained for _H. sapiens_ NA12878 gDNA, cDNA 1D2, and directRNA datasets, _Mus. musculus_ cDNA dataset, and the ZymoBIOMICS mock community datasets with 10 species and two abundance levels. Flowcell chemistry is R9.4 for all datasets. NA12878 gDNA and directRNA was basecalled by Guppy 3.1.5; NA12878 cDNA 1D2 was basecalled by Albacore 2.3.1; mouse cDNA was basecalled by Metrichor. These models are available within **[pre-trained_models folder](https://github.com/bcgsc/NanoSim/tree/master/pre-trained_models)**.  

### 2. Simulation stage
Simulation stage takes reference genome/transcriptome and read profiles as input and outputs simulated reads in FASTA format. Simulation stage runs in two modes: "genome" and "transcriptome" and you may use either of them based on your needs.

__Simulation stage general usage:__
```
usage: simulator.py [-h] [-v] {genome,transcriptome,metagenome} ...

Simulation step
-----------------------------------------------------------
Given error profiles, reference genome, metagenome,
and/or transcriptome, simulate ONT DNA or RNA reads

optional arguments:
  -h, --help            show this help message and exit
  -v, --version         show program's version number and exit

subcommands:
  
  There are two modes in read_analysis.
  For detailed usage of each mode:
      simulator.py mode -h
  -------------------------------------------------------

  {genome,transcriptome}
                        You may run the simulator on genome, transcriptome,
                        or metagenome mode.
    genome              Run the simulator on genome mode
    transcriptome       Run the simulator on transcriptome mode
    metagenome          Run the simulator on metagenome mode

```

**genome mode**  
If you are interested in simulating ONT genomic reads, you need to run the simulation stage in "genome" mode with following options.

__genome mode usage:__
```
usage: simulator.py genome [-h] -rg REF_G [-c MODEL_PREFIX] [-o OUTPUT]
                           [-n NUMBER] [-max MAX_LEN] [-min MIN_LEN]
                           [-med MEDIAN_LEN] [-sd SD_LEN] [--seed SEED]
                           [-k KMERBIAS] [-b {albacore,guppy,guppy-flipflop}]
                           [-s STRANDNESS] [-dna_type {linear,circular}]
                           [--perfect] [--fastq] [--chimeric] [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -rg REF_G, --ref_g REF_G
                        Input reference genome
  -c MODEL_PREFIX, --model_prefix MODEL_PREFIX
                        Location and prefix of error profiles generated from
                        characterization step (Default = training)
  -o OUTPUT, --output OUTPUT
                        Output location and prefix for simulated reads
                        (Default = simulated)
  -n NUMBER, --number NUMBER
                        Number of reads to be simulated (Default = 20000)
  -max MAX_LEN, --max_len MAX_LEN
                        The maximum length for simulated reads (Default =
                        Infinity)
  -min MIN_LEN, --min_len MIN_LEN
                        The minimum length for simulated reads (Default = 50)
  -med MEDIAN_LEN, --median_len MEDIAN_LEN
                        The median read length (Default = None), Note: this
                        simulation is not compatible with chimeric reads
                        simulation
  -sd SD_LEN, --sd_len SD_LEN
                        The standard deviation of read length in log scale
                        (Default = None), Note: this simulation is not
                        compatible with chimeric reads simulation
  --seed SEED           Manually seeds the pseudo-random number generator
  -k KMERBIAS, --KmerBias KMERBIAS
                        Minimum homopolymer length to simulate homopolymer
                        contraction and expansion events in, a typical k is 6
  -b {albacore,guppy,guppy-flipflop}, --basecaller {albacore,guppy,guppy-flipflop}
                        Simulate homopolymers with respect to chosen
                        basecaller: albacore, guppy, or guppy-flipflop
  -s STRANDNESS, --strandness STRANDNESS
                        Proportion of sense sequences. Overrides the value
                        profiled in characterization stage. Should be between
                        0 and 1
  -dna_type {linear,circular}
                        Specify the dna type: circular OR linear (Default =
                        linear)
  --perfect             Ignore error profiles and simulate perfect reads
  --fastq               Output fastq files instead of fasta files
  --chimeric            Simulate chimeric reads
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for simulation (Default = 1)

```

**transcriptome mode**  
If you are interested in simulating ONT transcriptome reads, you need to run the simulation stage in "transcriptome" mode with following options.

__transcriptome mode usage:__
```
usage: simulator.py transcriptome [-h] -rt REF_T [-rg REF_G] -e EXP
                                  [-c MODEL_PREFIX] [-o OUTPUT] [-n NUMBER]
                                  [-max MAX_LEN] [-min MIN_LEN] [--seed SEED]
                                  [-k KMERBIAS] [-b {albacore,guppy}]
                                  [-r {dRNA,cDNA_1D,cDNA_1D2}] [-s STRANDNESS]
                                  [--no_model_ir] [--perfect] [--polya POLYA]
                                  [--fastq] [-t NUM_THREADS] [--uracil]

optional arguments:
  -h, --help            show this help message and exit
  -rt REF_T, --ref_t REF_T
                        Input reference transcriptome
  -rg REF_G, --ref_g REF_G
                        Input reference genome, required if intron retention
                        simulation is on
  -e EXP, --exp EXP     Expression profile in the specified format as
                        described in README
  -c MODEL_PREFIX, --model_prefix MODEL_PREFIX
                        Location and prefix of error profiles generated from
                        characterization step (Default = training)
  -o OUTPUT, --output OUTPUT
                        Output location and prefix for simulated reads
                        (Default = simulated)
  -n NUMBER, --number NUMBER
                        Number of reads to be simulated (Default = 20000)
  -max MAX_LEN, --max_len MAX_LEN
                        The maximum length for simulated reads (Default =
                        Infinity)
  -min MIN_LEN, --min_len MIN_LEN
                        The minimum length for simulated reads (Default = 50)
  --seed SEED           Manually seeds the pseudo-random number generator
  -k KMERBIAS, --KmerBias KMERBIAS
                        Minimum homopolymer length to simulate homopolymer contraction and
                        expansion events in, a typical k is 6
  -b {albacore,guppy}, --basecaller {albacore,guppy}
                        Simulate homopolymers with respect to chosen  
                        basecaller: albacore or guppy
  -r {dRNA,cDNA_1D,cDNA_1D2}, --read_type {dRNA,cDNA_1D,cDNA_1D2}
                        Simulate homopolymers with respect to chosen read
                        type: dRNA, cDNA_1D or cDNA_1D2
  -s STRANDNESS, --strandness STRANDNESS
                        Proportion of sense sequences. Overrides the value
                        profiled in characterization stage. Should be between
                        0 and 1
  --no_model_ir         Ignore simulating intron retention events
  --perfect             Ignore profiles and simulate perfect reads
  --polya POLYA         Simulate polyA tails for given list of transcripts
  --fastq               Output fastq files instead of fasta files
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for simulation (Default = 1)
  --uracil              Converts the thymine (T) bases to uracil (U) in the
                        output fasta format
```

__sample expression file for transcriptome simulation__  
The expression profile is a tsv file containing expression levels of each isoform to be simulated. Users can use the output of `quantify` mode as template for modify or use the following format for constructing a new expression profile. 

The first row is header row specifying the format of the file. `target_id` is the id of each transcript. The `tpm` column is used for simulating, while the `est_count` is just a placeholder to be compatible with the output of the `quantify` mode and other quantification tools such as Salmon.

The following rows are entries for each transcript isoform and the id of which needs to exist in the provided reference transcriptome. The id should start with `ENS`.

Example:  
| target_id | est_counts | tpm |
| --------- |:----------:|----:|
| ENST00000222247.9 | 2307.2992 | 3145.3749 |  
| ENST00000274065.8 | 2641.9534 | 3601.5848 |  
| ENST00000400259.5 | 623.6130 | 850.1268 |  
| ENST00000344548.7 | 1828.3466 | 2492.4533 |  
| ENST00000484610.5 | 766.3528 | 1044.7137 |

**metagenome mode**  
If you are interested in simulating ONT metagenome reads, you need to run the simulation stage in "metagenome" mode with following options. We have provided sample config files for users to construct their own `-gl`, `-a`, and `-dl` config files correctly.

__metagenome mode usage:__
```
usage: simulator.py metagenome [-h] -gl GENOME_LIST -a ABUN -dl DNA_TYPE_LIST
                               [-c MODEL_PREFIX] [-o OUTPUT] [-max MAX_LEN]
                               [-min MIN_LEN] [-med MEDIAN_LEN] [-sd SD_LEN]
                               [--seed SEED] [-k KMERBIAS]
                               [-b {albacore,guppy,guppy-flipflop}]
                               [-s STRANDNESS] [--perfect]
                               [--abun_var ABUN_VAR [ABUN_VAR ...]] [--fastq]
                               [--chimeric] [-t NUM_THREADS]

optional arguments:
  -h, --help            show this help message and exit
  -gl GENOME_LIST, --genome_list GENOME_LIST
                        Reference metagenome list, tsv file, the first column
                        is species/strain name, the second column is the
                        reference genome fasta/fastq file directory
  -a ABUN, --abun ABUN  Abundance list, tsv file with header, the abundance of all species in
                        each sample need to sum up to 100. See example in README and provided
                        config files
  -dl DNA_TYPE_LIST, --dna_type_list DNA_TYPE_LIST
                        DNA type list, tsv file, the first column is
                        species/strain, the second column is the chromosome
                        name, the third column is the DNA type: circular OR
                        linear
  -c MODEL_PREFIX, --model_prefix MODEL_PREFIX
                        Location and prefix of error profiles generated from
                        characterization step (Default = training)
  -o OUTPUT, --output OUTPUT
                        Output location and prefix for simulated reads
                        (Default = simulated)
  -max MAX_LEN, --max_len MAX_LEN
                        The maximum length for simulated reads (Default =
                        Infinity)
  -min MIN_LEN, --min_len MIN_LEN
                        The minimum length for simulated reads (Default = 50)
  -med MEDIAN_LEN, --median_len MEDIAN_LEN
                        The median read length (Default = None), Note: this
                        simulationis not compatible with chimeric reads
                        simulation
  -sd SD_LEN, --sd_len SD_LEN
                        The standard deviation of read length in log scale
                        (Default = None), Note: this simulation is not
                        compatible with chimeric reads simulation
  --seed SEED           Manually seeds the pseudo-random number generator
  -k KMERBIAS, --KmerBias KMERBIAS
                        Minimum homopolymer length to simulate homopolymer
                        contraction andexpansion events in, a typical k is 6
  -b {albacore,guppy,guppy-flipflop}, --basecaller {albacore,guppy,guppy-flipflop}
                        Simulate homopolymers and/or base qualities with
                        respect to chosen basecaller: albacore, guppy, or
                        guppy-flipflop
  -s STRANDNESS, --strandness STRANDNESS
                        Percentage of antisense sequences. Overrides the value
                        profiled in characterization stage. Should be between
                        0 and 1
  --perfect             Ignore error profiles and simulate perfect reads
  --abun_var ABUN_VAR [ABUN_VAR ...]
                        Simulate random variation in abundance values, takes
                        in two values, format: relative_var_low,
                        relative_var_high, Example: -0.5 0.5)
  --fastq               Output fastq files instead of fasta files
  --chimeric            Simulate chimeric reads
  -t NUM_THREADS, --num_threads NUM_THREADS
                        Number of threads for simulation (Default = 1)
```

__sample abundance file for metagenome simulation__  

The abundance file is a tsv file, with rows representing the abundance of each sample and columns representing each sample. Each column (except for the first row) needs to sum up to 100, because the total abundance of each sample needs to be 100.  

The first row is header row to specify the number of reads in each sample. The format of the first row is:  
`Size  total_reads_in_sample1  total_reads_in_sample2  ...`  
The following rows are in the format as:  
`Species abundance_in_sample1  abundance_in_sample2  ...`    

| Size | 200000 | 100 |  
| ------------- |:-------------:| -----:|  
| Bacillus subtilis | 12 | 0.89 |  
| Cryptococcus neoformans | 2 | 0.00089 |  
| Enterococcus faecalis | 12 | 0.00089 |  
| Escherichia coli|	12 | 0.089 |  
| Lactobacillus fermentum |	12 |	0.0089 |  
| Listeria monocytogenes |	12 |	89.1 |  
| Pseudomonas aeruginosa	| 12 |	8.9 |  
| Saccharomyces cerevisiae |	2 |	0.89 |  
| Salmonella enterica	| 12 |	0.089 |  
| Staphylococcus aureus	| 12 |	0.000089 |  


In the above example, there are two samples. The first sample will contain 20,0000 reads, while the second sample will contain 100 reads. The abundances in sample 1 and 2 are as shown in th table, and both of them add up to 100.  

\* Notice: the use of `max_len` and `min_len` in genome mode will affect the read length distributions. If the range between `max_len` and `min_len` is too small, the program will run slowlier accordingly.  

\* Notice: the transcript name in the expression tsv file and the ones in th polyadenylated transcript list has to be consistent with the ones in the reference transcripts, otherwise the tool won't recognize them and don't know where to find them to extract reads for simulation.

\* Notice: the species name in the genome list file, dna type file, and abundance file has to be consistent. The chromosome names in the dna type file has to match the ones in the reference genomes. 

__Example runs:__  
1 If you want to simulate _E. coli_ genome, then circular command must be chosen because it's a circular genome  
`./simulator.py genome -dna_type circular -rg Ecoli_ref.fasta -c ecoli`

2 If you want to simulate only perfect reads, _i.e._ no snps, or indels, just simulate the read length distribution  
`./simulator.py genome -dna_type circular -rg Ecoli_ref.fasta -c ecoli --perfect`

3 If you want to simulate _S. cerevisiae_ genome with kmer bias, then linear command must be chosen because it's a linear genome  
`./simulator.py genome -dna_type linear -rg yeast_ref.fasta -c yeast --KmerBias`  

4 If you want to simulate human genome with length limits between 1000nt to 10000nt  
`./simulator.py genome -dna_type linear -rg human_ref.fasta -c human -max 10000 -min 1000`  

5 If you want to simulate human genome with controlled median read length and standard deviation, NanoSim will use log-normal distribution to approximate the read length distribution 
`./simulator.py genome -dna_type linear -rg human_ref.fasta -c human -med 5000 -sd 1.05`

6 If you want to simulate ten thousands cDNA/directRNA reads from mouse reference transcriptome  
`./simulator.py transcriptome -rt Mus_musculus.GRCm38.cdna.all.fa -rg Mus_musculus.GRCm38.dna.primary_assembly.fa -c mouse_cdna -e abundance.tsv -n 10000`

7 If you want to simulate five thousands cDNA/directRNA reads from mouse reference transcriptome without modeling intron retention  
`./simulator.py transcriptome -rt Mus_musculus.GRCm38.cdna.all.fa -c mouse_cdna -e abundance.tsv -n 5000 --no_model_ir`

8 If you want to simulate two thousands cDNA/directRNA reads from human reference transcriptome with polya tails, mimicking homopolymer bias (starting from homopolymer length >= 6) and reads in fastq format  
`./simulator.py transcriptome -rt Homo_sapiens.GRCh38.cdna.all.fa -c Homo_sapiens_model -e abundance.tsv -rg Homo_sapiens.GRCh38.dna.primary.assembly.fa --polya transcripts_with_polya_tails --fastq -k 6 --basecaller guppy -r dRNA`

9 If you want to simulate two metagenome samples with abundance variation, and chimeric reads  
`.simulator.py metagenome -gl sample_config_file/metagenome_list_for_simulation -a sample_config_file/abundance_for_simulation_multi_sample.tsv -dl sample_config_file/dna_type_list.tsv -c pre_trained_models/metagenome_ERR3152364_Even/training --abun_var -0.5 0.5 --chimeric`

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
15. `training_strandness_rate` Strandness rate in input reads

#### 1.2 Characterization stage (transcriptome)
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
16. `training_strandness_rate` Strandness rate in input reads
17. `training_addedintron_final.gff3` gff3 file format containing the intron coordinate information
18. `training_IR_info` List of transcripts in which there is a retained intron based on IR modeling step
19. `training_IR_markov_model` Markov model of Intron Retention events

#### 1.3 Characterization stage (metagenome)
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
15. `training_strandness_rate` Strandness rate in input reads
16. `training_chimeric_info` Information for chimeric reads  
17. `training_gap_length.pkl` The gap size between suplementary alignments  


### 2. Simulation stage  

1. `simulated_reads.fasta`
  FASTA file of simulated reads. Each reads has "unaligned", "aligned", or "perfect" in the header determining their error rate. "unaligned" means that the reads have an error rate over 90% and cannot be aligned. "aligned" reads have the same error rate as training reads. "perfect" reads have no errors.  
  
  To explain the information in the header, we have two examples:  
  * `>ref|NC-001137|-[chromosome=V]_468529_unaligned_0_F_0_3236_0`  
    All information before the first `_` are chromosome information. `468529` is the start position and `unaligned` suggesting it should be unaligned to the reference. The first `0` is the sequence index. `F` represents a forward strand. `0_3236_0` means that sequence length extracted from the reference is 3236 bases.  
  * `>ref|NC-001143|-[chromosome=XI]_115406_aligned_16565_R_92_12710_2`  
    This is an aligned read coming from chromosome XI at position 115406. `16565` is the sequence index. `R` represents a reverse complement strand. `92_12710_2` means that this read has 92-base head region (cannot be aligned), followed by 12710 bases of middle region, and then 2-base tail region.  
  
  The information in the header can help users to locate the read easily.  
  
__Specific to transcriptome simulation__: for reads that include retained introns, the header contains the information starting from `Retained_intron`, each genomic interval is separated by `;`.

__Specific to chimeric reads simulation__: for chimeric reads, different source chromosome and locations are separated by `;`, and there's a `chimeric` in the header to indicate.
  
2. `simulated_error_profile`
  Contains all the information of errors introduced into each reads, including error type, position, original bases and current bases.  


## Acknowledgements
Sincere thanks to our labmates and all contributors and users of this tool.
