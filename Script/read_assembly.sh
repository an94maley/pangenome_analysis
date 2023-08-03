## ./Script/read_assembly.sh
for each in `cat ./pinfo/Samples.txt`
do
seq_R1=./SPAdes/corrected_reads/${each}_R1.fastq.gz
seq_R2=./SPAdes/corrected_reads/${each}_R2.fastq.gz
spades.py -1 ${seq_R1} -2 ${seq_R2} -o ./SPAdes/assembly/${each}_result -t 100 -m 600 -k 65 --only-assembler --careful
mv ./SPAdes/assembly/${each}_result/scaffolds.fasta ./SPAdes/assembly/${each}_scaffolds.fasta
mv ./SPAdes/assembly/${each}_result/contigs.fasta ./SPAdes/assembly/${each}_contigs.fasta
done