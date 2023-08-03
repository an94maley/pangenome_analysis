## ./Script/rm_PCRdup.sh
for each in `cat ./pinfo/Samples.txt`
do
seq_R1=./Clean_data/${each}_R1.fastq
seq_R2=./Clean_data/${each}_R2.fastq
echo ${seq_R1} > ./Clean_data/${each}.list
echo ${seq_R2} >> ./Clean_data/${each}.list
fastuniq -i ./Clean_data/${each}.list -o ./Data_PCRrm/${each}_R1.fastq -p ./Data_PCRrm/${each}_R2.fastq
done