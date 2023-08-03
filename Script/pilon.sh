## ./Script/pilon.sh
Threads=30
mkdir BAM
mkdir Pilon
for each in `cat ./pinfo/Samples.txt`
do
path=./SPAdes/assembly/${each}_scaffolds.fasta
bwa index ${path}
bwa mem -t ${Threads} ${path} ./SPAdes/corrected_reads/${each}_R1.fastq.gz ./SPAdes/corrected_reads/${each}_R2.fastq.gz | samtools view -Sb - | samtools sort -@ ${Threads} -m 4G -o ./BAM/${each}.bam -O bam
samtools index -@ ${Threads} ./BAM/${each}.bam
pilon --genome ${path} --frags ./BAM/${each}.bam --output ${each} --outdir ./Pilon --threads ${Threads} --changes --vcf
done