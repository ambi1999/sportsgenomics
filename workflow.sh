#!/bin/bash


: <<'END'
##########commands for running wgsim
gcc -g -O2 -Wall -o wgsim wgsim.c -lz -lm
wgsim -N1000 -S1 RefGenome.fna simulated_reads/ambi_sim_reads.fq /dev/null
wgsim -N1000 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N7 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N7 -13 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N10 -13 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N8 -13 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N10 -13 -S1 RefGenome.fna ambi_sim_reads.fq /dev/null
wgsim -N10 -14 -S1 RefGenome.fna ambi_simulated_reads_len4.fq /dev/null
wgsim -N10 -14 -S1 AmbiExampleRefGenome.fna ambi_simulated_reads_len4.fq /dev/null

wgsim -N10 -14 -S1 Homo_sapiens_ACTN3_sequence.fa ambi_simulated_reads_ACTN3.fq /dev/null

wgsim -N1000 -S1 Homo_sapiens_ACTN3_sequence.fa ambi_simulated_reads_ACTN3.fq /dev/null

wgsim -N50000 -S1 Chr11GRCh38ACTN.fasta.fasta ambi_simulated_reads_ACTN3.fq /dev/null

/home/standalone/genomics/softwares/Sequencing read simulator/WGSIM/wgsim-master/Homo_sapiens_ACTN3_sequence_mutated.fa

wgsim -N20000 -S1 Homo_sapiens_ACTN3_sequence_mutated.fa ambi_simulated_reads_mutated_ACTN3.fq /dev/null

File: /media/standalone/OSDisk/ambi/Research/BOOK publication/Small Fastq and fasta Files/individual_1_ACTN3/ambi_simulated_1_reads_ACTN3_aligned.sorted.bam does not contain any sequence names which match the current genome.  File:      Ambi, Genome: chr1, chr2, chr3, chr4, ...

bowtie2-build AmbiExampleRefGenome1.fna  AmbiExampleRefGenome1

bowtie2-build AmbiExampleRefGenome1.fna  AmbiExampleRefGenome1

bowtie2-build Chr11GRCh38_ACN3.fasta Chr11GRCh38_ACN3

bowtie2-build Chr11GRCh38ACTN.fasta Chr11GRCh38ACTN

chr11.fa

bowtie2-build chr11.fa chr11.fa

Chr11GRCh38ACTN

/media/standalone/OSDisk/ambi/Research/BOOK publication/Small Fastq and fasta Files/ExampleGenome/AmbiExampleRefGenome.fna

/media/standalone/OSDisk/ambi/Research/BOOK publication/Small Fastq and fasta Files/Chr11Genome/Chr11GRCh38ACTN.fasta

bowtie2 -x ../ExampleGenome/AmbiExampleRefGenome -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam

bowtie2 -x ../Chr11Genome/Chr11GRCh38ACTN -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam

bowtie2 -x Chr11GRCh38/Chr11GRCh38ACTN -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam


END

: <<'END'

bowtie2-build ExampleGenome1/AmbiExampleRefGenome.fna  ExampleGenome1/AmbiExampleRefGenome
cd ExampleGenome1
ls -l
bowtie2-build AmbiExampleRefGenome1.fna  AmbiExampleRefGenome1
cd ..
bowtie2 -x ExampleGenome1/AmbiExampleRefGenome1 -U ambi_simulated_reads.fq -S alignments/ambi_simulated_reads_aligned.sam
ls -l
bowtie2 -x ExampleGenome1/AmbiExampleRefGenome1 -U ambi_simulated_reads.fq -S alignment/ambi_simulated_reads_aligned.sam
samtools view -b -S -o alignment/ambi_simulated_reads_aligned.sam alignment/ambi_simulated_reads_aligned.bam
samtools view -b -S -o alignment/ambi_simulated_reads_aligned.bam alignment/ambi_simulated_reads_aligned.sam
samtools sort alignment/ambi_simulated_reads_aligned.bam -o alignment/ambi_simulated_reads_aligned.sorted.bam
samtools index alignment/ambi_simulated_reads_aligned.sorted.bam
samtools view alignments/
samtools view alignment/ambi_simulated_reads_aligned.bam
samtools view alignment/ambi_simulated_reads_aligned.sorted.bam
samtools view alignment/ambi_simulated_reads_aligned.sam
history 1
history
history > commands.txt

END

: <<'END'


prefix="v1"

bowtie2 -x ../Chr11Genome/Chr11GRCh38ACTN -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam
samtools view -b -S -o ambi_simulated_1_reads_ACTN3_aligned.bam ambi_simulated_1_reads_ACTN3_aligned.sam
samtools sort ambi_simulated_1_reads_ACTN3_aligned.bam -o ambi_simulated_1_reads_ACTN3_aligned.sorted.bam
samtools index ambi_simulated_1_reads_ACTN3_aligned.sorted.bam
samtools view ambi_simulated_1_reads_ACTN3_aligned.sorted.bam

#bowtie2 -x /media/standalone/DATA/ambi/Human Genome/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam
END


: <<'END'

#bowtie2 -x "/media/standalone/DATA/ambi/HumanGenome/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index/GCA_000001405.15_GRCh38_full_analysis_set.fna.bowtie_index" -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam

#bowtie2 -x ../Chr11Genome/Chr11GRCh38ACTN -U ambi_simulated_1_reads_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam

bowtie2 -x ../Chr11Genome/chr11.fa -U ambi_simulated_reads_mutated_ACTN3.fq -S ambi_simulated_1_reads_ACTN3_aligned.sam
samtools view -b -S -o ambi_simulated_1_reads_ACTN3_aligned.bam ambi_simulated_1_reads_ACTN3_aligned.sam
samtools sort ambi_simulated_1_reads_ACTN3_aligned.bam -o ambi_simulated_1_reads_ACTN3_aligned.sorted.bam
samtools index ambi_simulated_1_reads_ACTN3_aligned.sorted.bam
#samtools view ambi_simulated_1_reads_ACTN3_aligned.sorted.bam

END

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
