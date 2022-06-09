
#!/bin/bash


source ~/anaconda3/etc/profile.d/conda.sh

conda activate pipits_env

pipits_funits
 
cd /home/nanoteam/Desktop/INTeGRATE_fastq_linux/ITS/Illumina/Fastq_Integrate_fungome
 
pispino_createreadpairslist -i /home/nanoteam/Desktop/INTeGRATE_fastq_linux/ITS/Illumina/Fastq_Integrate_fungome -o readpairslist.txt

pispino_seqprep -i /home/nanoteam/Desktop/INTeGRATE_fastq_linux/ITS/Illumina/Fastq_Integrate_fungome -o out_seqprep -l readpairslist.txt

pipits_funits -i out_seqprep/prepped.fasta -o out_funits -x ITS2

pipits_process -i out_funits/ITS.fasta -o out_process

pipits_funguild.py -i out_process/otu_table.txt -o out_process/otu_table_funguild.txt

history > his_09062022.txt
