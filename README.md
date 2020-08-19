# wga-pipeline

## New

### Pipeline

Launch this:

```
bash s01_wga_multisps.sh species.list genomes/ Tadh
```

* `species.list`: list of species to align (one per line)
* `genomes/`: folder containing all genomes to align (which have to have this format: Spsname_gDNA.fasta)
* `Tadh`: name of the reference species (all other species will be aligned to this one).

It will format genomes in the `genomes/` folder, split the reference genome into chromosomes (`Tadh`), index genomes, align genomes to reference genome, and post-process the alignments (see below for draft detailed pipeline).

### Dependencies

Alignets (`lastz` and maybe `multiz`?) and UCSC utilities can be installed from `bioconda`:

```bash
bioawk                    1.0                  hed695b0_5    bioconda
lastz                     1.0.4                h516909a_4    bioconda
multiz                    11.2                 h14c3975_0    bioconda
samtools                  1.10                 h9402c20_2    bioconda
ucsc-axtchain             377                  h446ed27_1    bioconda
ucsc-axttomaf             377                  h446ed27_1    bioconda
ucsc-bedsort              377                  h446ed27_2    bioconda
ucsc-chainnet             377                  h446ed27_1    bioconda
ucsc-chainprenet          377                  h446ed27_1    bioconda
ucsc-chainsort            377                  h446ed27_1    bioconda
ucsc-fatotwobit           377                  h446ed27_3    bioconda
ucsc-genepredtobed        377                  h446ed27_3    bioconda
ucsc-maftopsl             377                  h446ed27_1    bioconda
ucsc-netsyntenic          377                  h446ed27_1    bioconda
ucsc-nettoaxt             377                  h446ed27_1    bioconda
ucsc-nettobed             377                  h446ed27_1    bioconda
```

### Draft pipeline

Very rough draft pipeline for WGA alignment with **lastz** + **phastCons** + **phyloFit**.

*Note to XG: translate, clean, de-Anophelise, & annotate.*

1) Create genome indexes (fai and 2bit):

```
while read s ; do \
samtools faidx genomes_all/${s}_gDNA.fasta ;
/home/xavi/Programes/ucsc_utilities/faToTwoBit genomes_all/${s}_gDNA.fasta genomes_all/${s}_gDNA.fasta.2bit ;
done < ../00_mosq22sps.list
```

2) For the reference sps, split into chromosomes:

```
while read -a c ; do \
bioawk -c fastx '$1 == "'${c[0]}'" { print ">"$1"\n"$2 }' genomes_all/Anogam_gDNA.fasta \
| fold -w 80 > genomes_all/Anogam_gDNA.chr.${c[0]}.fasta ;
done < genomes_all/Anogam_gDNA.fasta.fai
```

3) MAF alignments (parameters as in Ensembl Compara; molt lent!):

	*https://www.ensembl.org/info/genome/compara/mlss.html?mlss=691*

```
mkdir alignments
while read -a c ; do \
while read s ; do \
lastz --step=10 --format=maf --gapped --gap=400,30 --seed=12of19 --transition --gappedthresh=3000 --inner=2200 --ambiguous=iupac --output=alignments/lastz_Anogam_${c[0]}-${s}.out.maf genomes_all/Anogam_gDNA.chr.${c[0]}.fasta genomes_all/${s}_gDNA.fasta \
done < 00_mosq22sps.list \
done < genomes_all/Anogam_gDNA.fasta.fai

```

4) Crear chains & nets a partir d'alineaments (lent!), i tornar a alineament MAF:

```
while read -a c ; do \
while read s ; do \
echo "##### ${c[0]}~${s} #####" ;
/home/xavi/Programes/ucsc_utilities/mafToPsl  $s Anogam alignments/lastz_Anogam_${c[0]}-${s}.out.maf alignments/lastz_Anogam_${c[0]}-${s}.out.psl ;
/home/xavi/Programes/ucsc_utilities/axtChain    -faQ -faT -linearGap=medium -psl alignments/lastz_Anogam_${c[0]}-${s}.out.psl genomes_all/Anogam_gDNA.fasta genomes_all/${s}_gDNA.fasta alignments/lastz_Anogam_${c[0]}-${s}.out.chain 2> /dev/null ;
/home/xavi/Programes/ucsc_utilities/chainSort   alignments/lastz_Anogam_${c[0]}-${s}.out.chain alignments/lastz_Anogam_${c[0]}-${s}.out.chain.sort ;
/home/xavi/Programes/ucsc_utilities/chainPreNet alignments/lastz_Anogam_${c[0]}-${s}.out.chain.sort genomes_all/Anogam_gDNA.fasta.fai genomes_all/${s}_gDNA.fasta.fai alignments/lastz_Anogam_${c[0]}-${s}.out.chain.filt ;
/home/xavi/Programes/ucsc_utilities/chainNet    alignments/lastz_Anogam_${c[0]}-${s}.out.chain.filt genomes_all/Anogam_gDNA.fasta.fai genomes_all/${s}_gDNA.fasta.fai stdout /dev/null | /home/xavi/Programes/ucsc_utilities/netSyntenic stdin alignments/lastz_Anogam_${c[0]}-${s}.out.net ;
/home/xavi/Programes/ucsc_utilities/netToBed    alignments/lastz_Anogam_${c[0]}-${s}.out.net alignments/lastz_Anogam_${c[0]}-${s}.out.net.bed ;
/home/xavi/Programes/ucsc_utilities/netToAxt    alignments/lastz_Anogam_${c[0]}-${s}.out.net alignments/lastz_Anogam_${c[0]}-${s}.out.chain.filt genomes_all/Anogam_gDNA.fasta.2bit genomes_all/${s}_gDNA.fasta.2bit alignments/lastz_Anogam_${c[0]}-${s}.out.net.axt ;
/home/xavi/Programes/ucsc_utilities/axtToMaf    alignments/lastz_Anogam_${c[0]}-${s}.out.net.axt genomes_all/Anogam_gDNA.fasta.fai genomes_all/${s}_gDNA.fasta.fai alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf ;
done < 00_mosq22sps.list ;
done < genomes_all/Anogam_gDNA.fasta.fai

```

	Steps:
	- mafToPsl, axtChain:   build chains from alignments
	- chainSort:            sort chains
	- chainNet,netSyntenic: build net from chains, filter syntenic nets
	- netToBed:             create simple bed file with netted regions in reference
	- netToAxt,axtToMaf:    create MAF output, good for phastCons
	- awk                   add sps name
	
	*Choices made*
	- in axtChain (chain creation), -linearGap is set to 'medium' instead of 'loose' 
	  (medium is better for closely related orgs)
	
	*RESULTING ALIGNMENTS DON'T HAVE OVERLAPS*
	
	

5a) Concatenate MAF files per reference chromosome & add sps reference

```
while read -a c ; do \
while read s ; do \
grep -v "^#" alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf \
| sed "1s/^/\n\n/" \
| awk 'BEGIN { RS="\n\na " ; FS="\ns " ; OFS=" " ; ORS="\n" } { print $1" "$2" "$3 }' \
| tr -s ' ' '\t' \
| awk 'NF > 0' \
| awk '{ print "a "$1"\ns Anogam."$2,$3,$4,$5,$6,$7"\ns '$s'."$8,$9,$10,$11,$12,$13"\n" }' ;
done < 00_mosq17sps.list \
| maf-sort > neti_Anogam_${c[0]}.maf ;
done < genomes_all/Anogam_gDNA.fasta.fai

```


*DEPRECATED*

```
while read -a c ; do \
while read s ; do \
grep -v "^#" alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf \
| sed "1s/^/\n\n/" \
| awk 'BEGIN { RS="\n\na " ; FS="\ns " ; OFS=" " ; ORS="\n" } { print $1" "$2" "$3 }' \
| tr -s ' ' '\t' \
| awk 'NF > 0' \
| awk '{ print $1,"Anogam."$2,$3,$4,$5,$6,$7,"'$s'."$8,$9,$10,$11,$12,$13 }' ;
done < 00_mosq22sps.list \
| sort -k 2,2 -k3,3n \
| awk '{ print "a "$1"\ns "$2,$3,$4,$5,$6,$7"\ns "$8,$9,$10,$11,$12,$13"\n" }' \
> nets_Anogam_${c[0]}.maf ;
done < genomes_all/Anogam_gDNA.fasta.fai

```


5b) Concatenate MAF files with MULTIZ:

```
while read -a c ; do \
while read s ; do \
if [ $s0 == $s ] ; then \
grep -v "^#" alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf \
| sed "1s/^/\n\n/" \
| awk 'BEGIN { RS="\n\na " ; FS="\ns " ; OFS=" " ; ORS="\n" } { print $1" "$2" "$3 }' \
| tr -s ' ' '\t' \
| awk 'NF > 0' \
| awk '{ print "a "$1"\ns Anogam."$2,$3,$4,$5,$6,$7"\ns '$s'."$8,$9,$10,$11,$12,$13"\n" }' \
| sed "1s/^/##maf version=1 scoring=blastz\n/" > mnet_Anogam_${c[0]}.out.maf
else grep -v "^#" alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf \
| sed "1s/^/\n\n/" \
| awk 'BEGIN { RS="\n\na " ; FS="\ns " ; OFS=" " ; ORS="\n" } { print $1" "$2" "$3 }' \
| tr -s ' ' '\t' \
| awk 'NF > 0' \
| awk '{ print "a "$1"\ns Anogam."$2,$3,$4,$5,$6,$7"\ns '$s'."$8,$9,$10,$11,$12,$13"\n" }' \
| sed "1s/^/##maf version=1 scoring=blastz\n/" > TMP ;
/home/xavi/Programes/multiz-tba.012109/multiz R=30 M=50 TMP mnet_Anogam_${c[0]}.out.maf 0 > mnet_Anogam_${c[0]}.out.maf2 && mv mnet_Anogam_${c[0]}.out.maf2 mnet_Anogam_${c[0]}.out.maf ;
fi ;
done < 00_mosq22sps_noQuery.list ;
maf-cut ${c[0]}:1-${c[1]} mnet_Anogam_${c[0]}.out.maf > mnet_Anogam_${c[0]}.out_clean.maf
done < genomes_all/Anogam_gDNA.fasta.fai

```

#### PhyloFit ####

6) Create training set with 100 alignments from each sps pair and each chromosome:

```
while read -a c ; do \
while read s ; do \
grep -v "^#" alignments/lastz_Anogam_${c[0]}-${s}.out.net.maf \
| sed "1s/^/\n\n/" \
| awk 'BEGIN { RS="\n\na " ; FS="\ns " ; OFS=" " ; ORS="\n" } { print $1" "$2" "$3 }' \
| tr -s ' ' '\t' \
| sed "s/score=/score /" \
| awk 'NF > 0 && $2 > 100' \
| sed "s/score /score=/" \
| shuf -n 100 \
| awk '{ print "a "$1"\ns Anogam."$2,$3,$4,$5,$6,$7"\ns '$s'."$8,$9,$10,$11,$12,$13"\n" }' ;
done < 00_mosq17sps.list \
| sed -r "s/s ([^_][^_][^_][^_][^_][^_])_[^ ]* /s \1 /" | maf-sort > phyloFit.chr.${c[0]}.training.maf
done < genomes_all/Anogam_gDNA.fasta.fai

```


7a) Phylofit to create general model (whole alignment):

```
mod=HKY ; \
while read -a c ; do \
phyloFit --tree 00_mosq17sps.newick --nrates 1 --subst-mod HKY85 --msa-format MAF --out-root phyloFit.chr.${c[0]}.${mod}.init phyloFit.chr.${c[0]}.training.maf
done < genomes_all/Anogam_gDNA.fasta.fai

```


#### PhastCons ####

8a) Use chromosome-specific models; initialize phastCons without priors (1) & launch estimate (2):

```
est=frei ; \
while read -a c ; do \
phastCons --estimate-trees phyloFit.chr.${c[0]}.${mod}.${est} --msa-format MAF phyloFit.chr.${c[0]}.training.maf phyloFit.chr.${c[0]}.${mod}.init.mod --no-post-probs ; \
phastCons --msa-format MAF neti_Anogam_${c[0]}.maf --viterbi neti_Anogam_${c[0]}.out.${mod}.${est}.viterbi.bed --score phyloFit.chr.${c[0]}.${mod}.${est}.cons.mod,phyloFit.chr.${c[0]}.${mod}.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${mod}.${est}.wig ;
sed "s/(null)/${c[0]}/" neti_Anogam_${c[0]}.out.${mod}.${est}.wig ;
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${mod}.${est}.wig ; \
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.cons.${mod}.${est}.wig neti_Anogam_00.out.cons.${mod}.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```


#### MOTIF ENRICHMENT: HOMER ####
*REVISE EXACT COMMANDS*


9) Extract intergenic CNEs (excluding those in exons and introns):

```
bedtools intersect -v -a neti_Anogam_00.out.scores.HKYG4.frei.viterbi.bed genomes_all/Anogam_long.annot.gff

```

10) Run HOMER to identify de novo motifs, and compare them to known motifs ()


```
findMotifsGenome.pl \
neti_Anogam_00.out.cons.HKYG4.frei.viterbi_nongenic.bed \
../genomes_all/Anogam_gDNA.fasta homer_consIntergenic_cisdbAg_CHNK_22ago18/ \
-mcheck /home/xavi/Programes/HOMER/motifs/DBs_manuals/Anogam_cisdb_20ago18.motif 
-size given -noknown -S 1000 \

```
Options:
* -noknown : doesn't check for known motifs, only denovo annotation (followed by simiilarity check)
* -mcheck  : motif database for simiilarity searches

Motifs from Drosophila gambiae Nitta et al 2015:
* /home/xavi/dades/Anotacions/tfbs_Nitta15_Elife/Dromel_Nitta15.homer.motif

Motifs from Anopheles gambiae in silico CIS-DB:
* /home/xavi/Programes/HOMER/motifs/DBs_manuals/Anogam_cisdb_20ago18.motif

Manual comparisons:


```
compareMotifs.pl homerMotifs.all.motifs comparama -known /home/xavi/dades/Anotacions/tfbs_CISBP_Weirauch14_20ago18/Anogam_cisdb_20ago18.motifs

```











#### Alternative steps ####



7b) For a *context-dependent method* you must specify less stringent algorithms:
```
mod=CDR2S ; \
while read -a c ; do \
phyloFit --tree 00_mosq17sps.newick --nrates 1 --subst-mod R2S --EM --precision MED --msa-format MAF --out-root phyloFit.chr.${c[0]}.${mod}.init phyloFit.chr.${c[0]}.training.maf
done < genomes_all/Anogam_gDNA.fasta.fai

```
8b) For a *context-dependent method* you must use SS alignments with tuples = N nucleotides in the model (2 or 3):
``` RUN ONCE:
while read -a c ; do \
msa_view phyloFit.chr.${c[0]}.training.maf --in-format MAF --tuple-size 2 --out-format SS > phyloFit.chr.${c[0]}.training.ss2
done < genomes_all/Anogam_gDNA.fasta.fai
```
```
est=frei ; \
while read -a c ; do \
phastCons --estimate-trees phyloFit.chr.${c[0]}.${mod}.${est} --msa-format SS phyloFit.chr.${c[0]}.training.ss3 phyloFit.chr.${c[0]}.${mod}.init.mod --no-post-probs ; \
phastCons --msa-format SS neti_Anogam_${c[0]}.ss3 --viterbi neti_Anogam_${c[0]}.out.${mod}.${est}.viterbi.bed --score phyloFit.chr.${c[0]}.${mod}.${est}.cons.mod,phyloFit.chr.${c[0]}.${mod}.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${mod}.${est}.wig ;
sed "s/chrom=[^ ]* /chrom=${c[0]} /" neti_Anogam_${c[0]}.out.${mod}.${est}.wig ;
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${mod}.${est}.wig ; \
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.${mod}.${est}.wig neti_Anogam_00.out.${mod}.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```





7c) PhyloFit to create neutral model from *4-fold degenerate sites* (vary freely):

```
while read -a c ; do \
awk '$1=="'${c[0]}'" && $3 =="CDS"' genomes_all/Anogam_long.annot.gff | sed "s/^/Anogam./" > genomes_all/Anogam_long.annot.chr.${c}.gff
msa_view phyloFit.chr.${c[0]}.training.maf   --in-format MAF --4d --features genomes_all/Anogam_long.annot.chr.${c[0]}.gff > phyloFit.chr.${c[0]}.training.4d.ss ; \
msa_view phyloFit.chr.${c[0]}.training.4d.ss --in-format SS  --out-format SS  --tuple-size 1 > phyloFit.chr.${c[0]}.training.4ds.ss ; \
phyloFit --tree 00_mosq17sps.newick --msa-format SS --out-root phyloFit.chr.${c[0]}.4dinit phyloFit.chr.${c[0]}.training.4ds.ss
done < genomes_all/Anogam_gDNA.fasta.fai

```
8c) Use chromosome-specific models with *4-fold degenerate sites*, estimating rho parameter and launching estimates in one run:
```
est=frei4dc ; \
while read -a c ; do \
phastCons --estimate-rho phyloFit.chr.${c[0]}.${est} --msa-format MAF phyloFit.chr.${c[0]}.training.maf phyloFit.chr.${c[0]}.4dinit.mod --no-post-probs ; \
phastCons --msa-format MAF neti_Anogam_${c[0]}.maf phyloFit.chr.${c[0]}.${est}.cons.mod,phyloFit.chr.${c[0]}.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${est}.wig ;
sed "s/(null)/${c[0]}/" neti_Anogam_${c[0]}.out.${est}.wig ; \
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${est}.wig ; \
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.${est}.wig neti_Anogam_00.out.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```






















8z) Use chromosome-specific models; initialize phastCons with priors (1) & launch estimate (2):

```
est=cov05len100 ; \
while read -a c ; do \
phastCons --target-coverage 0.05 --expected-length 100 --estimate-trees phyloFit.chr.${c[0]}.${est} --msa-format MAF phyloFit.chr.${c[0]}.training.maf phyloFit.chr.${c[0]}.init.mod --no-post-probs ; \
phastCons --msa-format MAF neti_Anogam_${c[0]}.maf phyloFit.chr.${c[0]}.${est}.cons.mod,phyloFit.chr.${c[0]}.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${est}.wig ;
sed "s/(null)/${c[0]}/" neti_Anogam_${c[0]}.out.${est}.wig ;
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${est}.wig ; \
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.${est}.wig neti_Anogam_00.out.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```





*PHASTCONS ALTERNATIVE: USE SPECIES-WIDE MODELS INSTEAD OF CHROMOSOME-SPECIFIC MODELS*

8x) Same sps-wide model; initialize phastCons without priors (1) & launch estimate (2):

```
est=frei ; \
phastCons --estimate-trees phyloFit.${est} --msa-format MAF phyloFit.training.maf phyloFit.v00.mod --no-post-probs ; \
while read -a c ; do \
phastCons --msa-format MAF neti_Anogam_${c[0]}.maf phyloFit.${est}.cons.mod,phyloFit.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${est}.wig ;
sed "s/(null)/${c[0]}/" neti_Anogam_${c[0]}.out.${est}.wig ;
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${est}.wig ;
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.${est}.wig neti_Anogam_00.out.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```

8y) Same sps-wide model; Initialize phastCons with priors (1) & launch estimate (2):

```
est=cov70len20 ; \
phastCons --target-coverage 0.70 --expected-length 20 --estimate-trees phyloFit.${est} --msa-format MAF phyloFit.training.maf phyloFit.v00.mod --no-post-probs ; \
while read -a c ; do \
phastCons --msa-format MAF neti_Anogam_${c[0]}.maf phyloFit.${est}.cons.mod,phyloFit.${est}.noncons.mod > neti_Anogam_${c[0]}.out.${est}.wig ;
sed "s/(null)/${c[0]}/" neti_Anogam_${c[0]}.out.${est}.wig ;
done < genomes_all/Anogam_gDNA.fasta.fai > neti_Anogam_00.out.${est}.wig ;
/home/xavi/Programes/IGVTools/igvtools toTDF neti_Anogam_00.out.${est}.wig neti_Anogam_00.out.${est}.tdf genomes_all/Anogam_gDNA.fasta.fai

```
