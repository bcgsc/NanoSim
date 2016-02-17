#!/bin/bash

#####
# This script shows how the simulated reads are generated in the paper.
# All training datasets are downloaded from ENA and processed with poretools.
# Only 2D pass reads are extracted and stored on bcgsc ftp server.
#####

# Download the source file
mkdir NanoSim
cd NanoSim
wget https://github.com/bcgsc/NanoSim/archive/master.zip


# Create a working directory
mkdir ecoli_simulation
cd ecoli_simulation

# 1. E. coli R7 dataset
# Origin: ftp://climb.genomics.cn/pub/10.5524/100001_101000/100102/Ecoli_R7_CombinedFasta.tgz

# Get the 2D reads
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/*
# Get the reference genome
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/*

# Profiling stage


# To get the profile directly:
wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/ecoli_R7_profile.zip

