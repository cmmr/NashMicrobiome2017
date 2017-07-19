#!/bin/sh

THREADS=$1;

if [ -z "$THREADS" ];
	then THREADS=`grep -c ^processor /proc/cpuinfo`;
fi



cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_mergepairs {}.1.fq -reverse {}.2.fq -fastqout {}.merged.fq -fastq_allowmergestagger -fastq_minovlen 50';
cat ../SampleList | parallel -I {} 'echo {}; perl ~mcwong/mergeReads.pl {}.1.fq {}.2.fq {}.merged.fq > {}.temp.fq';
cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; usearch70 -fastq_filter {}.temp.fq -relabel "{}_" -fastqout {}.filtered.fq -eeout -fastq_minlen 200';
cat ../SampleList | parallel -j$THREADS -I {} 'echo {}; perl ~mcwong/filterSeqs.pl {}.filtered.fq > {}.filtered2.fq'
cat *.filtered2.fq > seqs.fq;
bowtie2 -x /gpfs1/db/phix/bowtie2/phix  -U seqs.fq --end-to-end --very-sensitive -p $THREADS --un seqs.filtered.fq -S /dev/null 2>../phix.bleed.txt;
mkdir ../split_libraries;
fq2fa seqs.filtered.fq seqs.fna;
mv seqs.fna ../split_libraries;
rm *.filtered2.fq;
rm *.temp.fq;
rm *.merged.fq;
rm *.filtered.fq;
rm seqs.fq;
rm seqs.filtered.fq;