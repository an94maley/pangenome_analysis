## ./Script/QC_trim.sh
Trimmomatic="/d2/home/user3/miniconda3/ENTER/pkgs/trimmomatic-0.36-3/share/trimmomatic-0.36-3/trimmomatic.jar"
SampleID=$1
THREADS=30
mkdir Log
mkdir Clean_data
mkdir -p Trim/${SampleID}_trimming
cd Trim/${SampleID}_trimming
ln -s ../../Raw_data/
ln -s ../../Log/
ln -s ../../Clean_data/
SM=${SampleID}
LB=${SampleID}
ID=${SampleID}
PU=${SampleID}
RGinfo="@RG\tID:${ID}\tSM:${SM}\tPL:illumina\tLB:${LB}\tPU:${PU}"
f1=./Raw_data/${SampleID}_R1.fastq.gz
f2=./Raw_data/${SampleID}_R2.fastq.gz
java -jar ${Trimmomatic} PE -threads ${THREADS} ${f1} ${f2} R1_paired_fastq.gz R1_unpaired_fastq.gz R2_paired_fastq.gz R2_unpaired_fastq.gz ILLUMINACLIP:/d2/home/user3/miniconda3/ENTER/pkgs/trimmomatic-0.36-3/share/trimmomatic-0.36-3/adapters/TruSeq3-PE.fa:2:30:10 LEADING:5 TRAILING:5 SLIDINGWINDOW:5:20 HEADCROP:10 MINLEN:75 TOPHRED33 2> ./Log/${SampleID}_QC_clean.log
cp R1_paired_fastq.gz ./Clean_data/${SampleID}_R1.fastq.gz
cp R2_paired_fastq.gz ./Clean_data/${SampleID}_R2.fastq.gz
rm -r ../${SampleID}_trimming
cd ../..