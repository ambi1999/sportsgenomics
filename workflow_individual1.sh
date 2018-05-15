#!/bin/bash


: <<'END'

echo "Generating artificial reads"
wgsim -N20000 -S1 Homo_sapiens_ACTN3_sequence_mutated.fa ambi_simulated_reads_mutated_ACTN3.fq /dev/null


File: /media/standalone/OSDisk/ambi/Research/BOOK publication/Small Fastq and fasta Files/individual_1_ACTN3/ambi_simulated_1_reads_ACTN3_aligned.sorted.bam does not contain any sequence names which match the current genome.  File:      Ambi, Genome: chr1, chr2, chr3, chr4, ...

bowtie2-build AmbiExampleRefGenome1.fna  AmbiExampleRefGenome1

bowtie2-build AmbiExampleRefGenome1.fna  AmbiExampleRefGenome1

bowtie2-build Chr11GRCh38_ACN3.fasta Chr11GRCh38_ACN3

bowtie2-build Chr11GRCh38ACTN.fasta Chr11GRCh38ACTN

chr11.fa

echo "Build index for alignment"

bowtie2-build chr11.fa chr11

END


: <<'END'

bowtie2 -x ../Chr11Genome/chr11 -U individual_1_reads.fq -S individual_1_aligned.sam
samtools view -b -S -o individual_1_aligned.bam individual_1_aligned.sam
samtools sort individual_1_aligned.bam -o individual_1_aligned.sorted.bam
samtools index individual_1_aligned.sorted.bam
samtools view individual_1_aligned.sorted.bam



prefix="mutated_ACTN3"
#pathVarscan="/home/ubuntu/ambi/software/varscan/"
#/home/standalone/genomics/softwares/Varscan/VarScan.v2.4.0.jar
pathVarscan="/home/standalone/genomics/softwares/Varscan/"

#echo "generate vcf using freebayes"\
#&& freebayes -f ../Chr11Genome/chr11.fa ambi_simulated_1_reads_ACTN3_aligned.sorted.bam > "$prefix"_variants_freebayes.vcf\
echo "generate mpileup"\
&&samtools mpileup -f ../Chr11Genome/chr11.fa ambi_simulated_1_reads_ACTN3_aligned.sorted.bam > "$prefix".mpileup\
&&echo "generate vcf using varscan and mpileup using varscan mpileup2cns"\
&&java -jar $pathVarscan/VarScan.v2.4.0.jar mpileup2cns "$prefix".mpileup --min-coverage 40 --min-reads2 20 --output-vcf 1 --variants > "$prefix"_variants_varscan.vcf\

END

pathVarscan="/home/standalone/genomics/softwares/Varscan/"

bowtie2 -x ../Chr11Genome/chr11 -U individual_1_reads.fq -S individual_1_aligned.sam

samtools view -b -S -o individual_1_aligned.bam individual_1_aligned.sam

samtools sort individual_1_aligned.bam -o individual_1_aligned.sorted.bam

samtools index individual_1_aligned.sorted.bam

samtools view individual_1_aligned.sorted.bam

samtools mpileup -f ../Chr11Genome/chr11 individual_1_aligned.sorted.bam > individual_1_.mpileup\

pathVarscan="/home/standalone/genomics/softwares/Varscan/"

java -jar $pathVarscan/VarScan.v2.4.0.jar mpileup2cns individual_1_.mpileup --min-coverage 40 --min-reads2 20 --output-vcf 1 --variants > individual_1_variants_varscan.vcf

freebayes -f ../Chr11Genome/chr11.fa individual_1_aligned.sorted.bam > individual_1_variants_freebayes.vcf\

prefix="mutated_ACTN3"
#pathVarscan="/home/ubuntu/ambi/software/varscan/"
#/home/standalone/genomics/softwares/Varscan/VarScan.v2.4.0.jar
pathVarscan="/home/standalone/genomics/softwares/Varscan/"

#echo "generate vcf using freebayes"\
#&& freebayes -f ../Chr11Genome/chr11.fa ambi_simulated_1_reads_ACTN3_aligned.sorted.bam > "$prefix"_variants_freebayes.vcf\
echo "generate mpileup"\
&&samtools mpileup -f ../Chr11Genome/chr11.fa ambi_simulated_1_reads_ACTN3_aligned.sorted.bam > "$prefix".mpileup\
&&echo "generate vcf using varscan and mpileup using varscan mpileup2cns"\
&&java -jar $pathVarscan/VarScan.v2.4.0.jar mpileup2cns "$prefix".mpileup --min-coverage 40 --min-reads2 20 --output-vcf 1 --variants > "$prefix"_variants_varscan.vcf\



