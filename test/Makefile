.PHONY: all clean

all: ecoli_R7_2D.fasta ecoli_K12_MG1655_ref.fa
	nanosimh_train -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -o model_ecoli -b 1
	nanosimh_simulate -r ecoli_K12_MG1655_ref.fa -c model_ecoli -n 1 --circular
	

ecoli_R7_2D.fasta:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 200 > $@

ecoli_K12_MG1655_ref.fa:
	wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@


clean:
	rm -f *.fa *.maf
	
