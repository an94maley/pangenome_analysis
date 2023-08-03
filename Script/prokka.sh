## ./Script/prokka.sh
mkdir Prokka
for each in `cat ./pinfo/Samples.txt`
do
path=./Pilon/${each}.fasta
prokka --genus Mycobacterium --kingdom Bacteria --cpus 16 --outdir ./Prokka/${each} --prefix ${each} --compliant --force ${path}
done