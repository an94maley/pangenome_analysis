## ./Script/read_correct.sh
for each in `cat ./pinfo/Samples.txt`
do
seq_R1=./Data_PCRrm/${each}_R1.fastq.gz
seq_R2=./Data_PCRrm/${each}_R2.fastq.gz
spades.py -1 ${seq_R1} -2 ${seq_R2} -o ./SPAdes/corrected_reads -t 100 -m 600 --only-error-correction
mv ./SPAdes/corrected_reads/corrected/${each}_R1.fastq.00.0_0.cor.fastq.gz ./SPAdes/corrected_reads/${each}_R1.fastq.gz
mv ./SPAdes/corrected_reads/corrected/${each}_R2.fastq.00.0_0.cor.fastq.gz ./SPAdes/corrected_reads/${each}_R2.fastq.gz
done