.PHONY: all clean

srcdir=../src

all: ecoli_R7_2D.fasta ecoli_K12_MG1655_ref.fa
	$(srcdir)/nanosim-train.py -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -o ecoli -b 1
	$(srcdir)/nanosim-simulate.py circular -r ecoli_K12_MG1655_ref.fa -c ecoli -n 1 # Note the -c option has to be the same as -o in read_analysis.py, or both use default parameter
	

ecoli_R7_2D.fasta:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 100 > $@

ecoli_K12_MG1655_ref.fa:
	wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@


clean:
	rm -f *.fa *.maf
	