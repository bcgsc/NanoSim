.PHONY: all clean train simulate

all: train simulate

train: model_ecoli_model_profile

model_ecoli_model_profile: ecoli_R7_2D.fasta ecoli_K12_MG1655_ref.fa
	nanosimh_train -i ecoli_R7_2D.fasta -r ecoli_K12_MG1655_ref.fa -p model_ecoli

simulate: ecoli_K12_MG1655_ref.fa
	nanosimh_simulate -r ecoli_K12_MG1655_ref.fa -p model_ecoli -n 1 --circular -o sim1
	nanosimh_simulate -r ecoli_K12_MG1655_ref.fa -p model_ecoli -n 1 --circular -o sim2
	diff sim1_reads.fasta sim2_reads.fasta

ecoli_R7_2D.fasta:
	curl ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@ | head -n 200 > $@

ecoli_K12_MG1655_ref.fa:
	wget ftp://ftp.bcgsc.ca/supplementary/NanoSim/$@

clean:
	rm -f *.fa *.maf model_* ref_genome.* *.log *.fasta simulated_error_profile

