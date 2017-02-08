#!/bin/bash

#####
# This script shows how the simulated reads are generated in the paper.
# All training datasets are downloaded from ENA and processed with poretools.
# Only 2D pass reads are extracted and stored on bcgsc ftp server.
#####

##### Download the source file
mkdir NanoSim
cd NanoSim
wget https://github.com/bcgsc/NanoSim/archive/master.zip
unzip master.zip
# After this step, you should have a folder called NanoSim, and inside you have master.zip and NanoSim-master two sub-folders

##### Inside NanoSim, create a working directory
mkdir ecoli_simulation
cd ecoli_simulation

# 1. E. coli R7 dataset
# Origin: ftp://climb.genomics.cn/pub/10.5524/100001_101000/100102/Ecoli_R7_CombinedFasta.tgz

# Get the 2D reads
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_R7_2D.fasta
# Get the reference genome
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_K12_MG1655_ref.fa

# Profiling stage, make sure to set the mode of read_analysis.py to -r-x or above
../NanoSim-master/src/read_analysis.py -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -o ecoli

# Simulation stage, suppose the genome to be simulated is called test.fasta and make sure to provide the correct path to it
../NanoSim-master/src/simulator.py circular -r test.fasta -c ecoli # Note the -c option has to be the same as -o in read_analysis.py, or both use default parameter

# To get the profile directly:
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_R7_profile.zip

# 2. E. coli R7.3 dataset
# Origin: http://www.ebi.ac.uk/ena/data/view/ERX708228, ERX708229, ERX708230, ERX708231

# Get the 2D reads
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_R73_2D.fasta
# Get the reference genome
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_K12_MG1655_ref.fa

# To get the profile directly:
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_R73_profile.zip

# 3. E. coli UCSC phase1b dataset
# Origin:  http://www.ebi.ac.uk/ena/data/view/ERP010368

# Get the 2D reads
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_UCSC_phase1b_2D.fasta
# Get the reference genome
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_K12_MG1655_ref.fa

# To get the profile directly:
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_UCSC1b_profile.zip

# 4. S. cerevisiae dataset
# Origin: http://labshare.cshl.edu/shares/schatzlab/www-data/nanocorr/2015.07.07/W303_ONT_Raw_reads.fa.gz

# Get the 2D reads
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/yeast_2D.fasta
# Get the reference genome
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/yeast_S288C_ref.fa

# To get the profile directly:
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/yeast_profile.zip

