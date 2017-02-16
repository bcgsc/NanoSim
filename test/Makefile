.PHONY: all clean

srcdir=../src


all: ref.fa reads.fa

	$(srcdir)/read_analysis.py -i reads.fa -r ref.fa -o ecoli -b 1
	$(srcdir)/simulator.py circular -r ref.fa -c ecoli -n 1

#all: ecoli_R7_2D.fasta ecoli_K12_MG1655_ref.fa
#	$(srcdir)/read_analysis.py -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -o ecoli -b 1
#	$(srcdir)/simulator.py circular -r ecoli_K12_MG1655_ref.fa -c ecoli -n 1 # Note the -c option has to be the same as -o in read_analysis.py, or both use default parameter
	

ecoli_R7_2D.fasta:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 100 > $@

ecoli_K12_MG1655_ref.fa:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 3 > $@


clean:
	rm -f *.fa *.maf
	
