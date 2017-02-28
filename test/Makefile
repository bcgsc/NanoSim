.PHONY: all clean

all: ecoli_R7_2D.fasta ecoli_K12_MG1655_ref.fa
	nanosimh_train -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -p model_ecoli
	nanosimh_simulate -r ecoli_K12_MG1655_ref.fa -p model_ecoli -n 1 --circular
	

ecoli_R7_2D.fasta:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 200 > $@

ecoli_K12_MG1655_ref.fa:
	wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@


clean:
	rm -f *.fa *.maf model_* ref_genome.* *.log *.fasta
