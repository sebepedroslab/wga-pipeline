#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q long-sl7
#$ -l virtual_free=10G,h_rt=720:00:00
#$ -o tmp/
#$ -e tmp/

spslist=$1
gen_fn=$2
ref=$3
ref_min_length=500000
conda_path=$(which lastz)
conda_path=$(dirname $conda_path)

# polish input genomes
mkdir -p query_genomes/
while read i ; do

	if [ -s query_genomes/${i}_gDNA.fasta.2bit ]
		then echo "already clean ${i}"
		else echo "clean input ${i}"
			sed "s/>/>${i}./" ${gen_fn}/${i}_gDNA.fasta > query_genomes/${i}_gDNA.fasta
			samtools faidx query_genomes/${i}_gDNA.fasta
			faToTwoBit query_genomes/${i}_gDNA.fasta query_genomes/${i}_gDNA.fasta.2bit
	fi

done < ${spslist}

# split reference into chromosomes
mkdir -p ref_chromosomes/
bioawk -c fastx 'length($2) > '"${ref_min_length}"' { print ">"$1"\n"$2 }' query_genomes/${ref}_gDNA.fasta > ref_${ref}_only_long_chr.fasta
bioawk -c fastx 'length($2) > '"${ref_min_length}"' { print $1 }'          query_genomes/${ref}_gDNA.fasta > ref_${ref}_only_long_chr.list.txt
while read c ; do
	if [ -s ref_chromosomes/ref_${ref}_gDNA.${c}.fasta ]
		then echo "already split: ${ref} chromosome ${c[0]}"
		else echo "split: ${ref} chromosome ${c[0]}"
			bioawk -c fastx '$1 == "'${c}'" { print ">"$1"\n"$2 }' ref_${ref}_only_long_chr.fasta | fold -w 200 > ref_chromosomes/ref_${ref}_gDNA.${c}.fasta
			samtools faidx ref_chromosomes/ref_${ref}_gDNA.${c}.fasta
			faToTwoBit ref_chromosomes/ref_${ref}_gDNA.${c}.fasta ref_chromosomes/ref_${ref}_gDNA.${c}.fasta.2bit
	fi
done < ref_${ref}_only_long_chr.list.txt

# alignments to reference species

# function: lastz alignments and downstream parsing
do_lastz="
#!/bin/bash
#$ -V
#$ -cwd
#$ -M xavier.graubove@crg.eu
#$ -m a
#$ -q mem_512_12h
#$ -l virtual_free=50G,h_rt=43200
#$ -o tmp/
#$ -e tmp/

refsps=\$1
reffas=\$2
qrysps=\$3
qryfas=\$4
output=\$5

# do alignments
${conda_path}/lastz \
--step=10 --format=maf --gapped --gap=400,30 --seed=12of19 \
--transition --gappedthresh=3000 --inner=2200 --ambiguous=iupac \
--output=\${output}.maf \
\${reffas} \
\${qryfas}

# alter alignment to ensure seq names are conserved
sed -i \"s/^s \${refsps}\./s \${refsps}.\${refsps}./\" \${output}.maf
sed -i \"s/^s \${qrysps}\./s \${qrysps}.\${qrysps}./\" \${output}.maf

# polish alignments
${conda_path}/mafToPsl    \${qrysps} \${refsps} \${output}.maf \${output}.psl
${conda_path}/axtChain    -faQ -faT -linearGap=medium -psl \${output}.psl \${reffas} \${qryfas} \${output}.chain
${conda_path}/chainSort   \${output}.chain \${output}.chain.sort
${conda_path}/chainPreNet \${output}.chain.sort \${reffas}.fai \${qryfas}.fai \${output}.chain.filt
${conda_path}/chainNet    \${output}.chain.filt \${reffas}.fai \${qryfas}.fai stdout /dev/null | ${conda_path}/netSyntenic stdin \${output}.net
${conda_path}/netToBed    \${output}.net \${output}.net.bed
${conda_path}/netToAxt    \${output}.net \${output}.chain.filt \${reffas}.2bit \${qryfas}.2bit \${output}.net.axt
${conda_path}/axtToMaf    \${output}.net.axt \${reffas}.fai \${qryfas}.fai \${output}.net.maf
"


mkdir -p tmp/
mkdir -p maf_alignments/
while read c ; do
while read i ; do

	if [[ $i != ${ref} ]] ; then
		echo "align: ${i} to ${ref} ${c}"
		echo "${do_lastz}" > tmp/qsub_annotation.${ref}.${c}-${i}.sh
		qsub -N lastz.${ref}.${c}-${i} tmp/qsub_annotation.${ref}.${c}-${i}.sh \
			${ref} \
			ref_chromosomes/ref_${ref}_gDNA.${c}.fasta \
			${i}  \
			query_genomes/${i}_gDNA.fasta \
			maf_alignments/lastz_${ref}.${c}-${i}.out
	fi

done < ${spslist}
done < ref_${ref}_only_long_chr.list.txt




exit

while read -a c ; do
while read i ; do

	if [[ $i != ${ref} ]] ; then
		echo "process alignments: ${i} to ${ref} ${c}"
		mafToPsl    ${i} ${ref} maf_alignments/lastz_${ref}.${c}-${i}.out.maf maf_alignments/lastz_${ref}.${c}-${i}.out.psl
		axtChain    -faQ -faT -linearGap=medium -psl maf_alignments/lastz_${ref}.${c}-${i}.out.psl ${gen_fn}/${ref}_gDNA.fasta ${gen_fn}/${i}_gDNA.fasta maf_alignments/lastz_${ref}.${c}-${i}.out.chain
		chainSort   maf_alignments/lastz_${ref}.${c}-${i}.out.chain maf_alignments/lastz_${ref}.${c}-${i}.out.chain.sort
		chainPreNet maf_alignments/lastz_${ref}.${c}-${i}.out.chain.sort ${gen_fn}/${ref}_gDNA.fasta.fai ${gen_fn}/${i}_gDNA.fasta.fai maf_alignments/lastz_${ref}.${c}-${i}.out.chain.filt
		chainNet    maf_alignments/lastz_${ref}.${c}-${i}.out.chain.filt ${gen_fn}/${ref}_gDNA.fasta.fai ${gen_fn}/${i}_gDNA.fasta.fai stdout /dev/null | netSyntenic stdin maf_alignments/lastz_${ref}.${c}-${i}.out.net
		netToBed    maf_alignments/lastz_${ref}.${c}-${i}.out.net maf_alignments/lastz_${ref}.${c}-${i}.out.net.bed
		netToAxt    maf_alignments/lastz_${ref}.${c}-${i}.out.net maf_alignments/lastz_${ref}.${c}-${i}.out.chain.filt ${gen_fn}/${ref}_gDNA.fasta.2bit ${gen_fn}/${i}_gDNA.fasta.2bit maf_alignments/lastz_${ref}.${c}-${i}.out.net.axt
		axtToMaf    maf_alignments/lastz_${ref}.${c}-${i}.out.net.axt ${gen_fn}/${ref}_gDNA.fasta.fai ${gen_fn}/${i}_gDNA.fasta.fai maf_alignments/lastz_${ref}.${c}-${i}.out.net.maf
	fi

done < ${spslist}
done < ref_${ref}_only_long_chr.list.txt
