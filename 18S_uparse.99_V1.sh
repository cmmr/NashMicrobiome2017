#!/bin/sh

THREADS=$1;

if [ -z "$THREADS" ];
        then THREADS=`grep -c ^processor /proc/cpuinfo`;
fi


mkdir uparse;
usearch70 -derep_fulllength split_libraries/seqs.fna  -output uparse/derep.fna -sizeout -uc uparse/derep.uc -threads $THREADS;
usearch70 -sortbysize       uparse/derep.fna -output uparse/sorted.fa -minsize 2;
cp uparse/sorted.fa uparse/temp.fa
for i in {0.4,0.7,1.0};
do 
usearch70 -cluster_otus     uparse/temp.fa -otus   uparse/temp1.fa -otu_radius_pct $i -uc uparse/cluster_$i.uc -fastaout uparse/clustering.$i.fasta.out;
cat uparse/clustering.$i.fasta.out | grep "^>" | grep chimera | sed 's/^>//g' | sed -re 's/;n=.*up=/\t/g' | sed 's/;$//g' | tee -a uparse/chimeras.txt > uparse/chimeras.$i.txt;
cat uparse/clustering.$i.fasta.out | grep "^>" > uparse/uparseref.decisions.$i.txt;
rm uparse/clustering.$i.fasta.out;
mv uparse/temp1.fa uparse/temp.fa;
done
mv uparse/temp.fa uparse/otus1.fa
usearch80 -uchime_ref       uparse/otus1.fa  -db /gpfs1/projects/auchtung/silva123.usearch80.udb -strand plus -uchimeout uparse/uchimeref.uc -threads $THREADS
cat uparse/uchimeref.uc | cut -f2,17 | grep -v "Y$" | cut -f1 | ~mcwong/getSeq uparse/otus1.fa > uparse/otus.fa;
usearch80 -usearch_global uparse/otus.fa -db /gpfs1/projects/auchtung/silva123.usearch80.udb -maxaccepts 0 -maxrejects 0 -strand both -id .97 -query_cov .95 -threads $THREADS -uc uparse/97.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
perl ~mcwong/GitRepos/mcwong/ITS_workflows/cleanHitsTableITS.pl uparse/97.centroids.uc /gpfs1/db/silva/123/silva123.map > uparse/97.clean.uc;
mv NewTaxa.txt uparse/NewTaxa.97.txt;
cat uparse/97.centroids.uc | grep "*$" | cut -f9 | ~mcwong/getSeq uparse/otus.fa > uparse/missed.fa;
usearch80 -usearch_global uparse/missed.fa -db /gpfs1/projects/auchtung/silva123.usearch80.udb -maxaccepts 0 -maxrejects 0 -strand both -id .80 -query_cov .95 -threads $THREADS -uc uparse/missed.centroids.uc -gapopen 5.0I/0.0E -gapext 1.0I/0.0E -top_hits_only;
cat uparse/97.clean.uc | grep -v "*$" | sed -re 's/\t1$//g' | cut -f10 | sed 's/.*/*/g' > second;
cat uparse/97.clean.uc | grep -v "*$" | sed -re 's/\t1$//g' | cut -f1-9 > first;
paste first second >> uparse/missed.centroids.uc;
rm first second
perl ~mcwong/GitRepos/mcwong/ITS_workflows/cleanHitsTableITS.pl uparse/missed.centroids.uc /gpfs1/db/silva/123/silva123.map > uparse/missed.clean.uc;
mv NewTaxa.txt uparse/NewTaxa.missed.txt
cat uparse/derep.fna | grep -A1 "size=1;" | cut -f2 -d ">" | ~mcwong/getSeq uparse/derep.fna > uparse/singletons.fna;
usearch70 -usearch_global uparse/singletons.fna -db uparse/sorted.fa -id .99 -uc uparse/singletons2otus.uc -strand plus -threads $THREADS -maxaccepts 32 -maxrejects 128 -query_cov .85 -wordlength 12;
perl ~mcwong/GitRepos/mcwong/16S_workflows/resolveIterativeUparse.pl uparse/cluster_* uparse/singletons2otus.uc uparse/97.clean.uc --derep uparse/derep.uc --chimeras uparse/chimeras.txt --taxonomy uparse/NewTaxa.97.txt --uchime uparse/uchimeref.uc
mv otu_table.biom uparse/otu_table_97.biom;
mv reads2otus.txt uparse/reads2otus.97.txt;
perl ~mcwong/GitRepos/mcwong/16S_workflows/resolveIterativeUparse.pl uparse/cluster_* uparse/singletons2otus.uc uparse/missed.clean.uc --perc --otus uparse/missed.clean.uc  --derep uparse/derep.uc --chimeras uparse/chimeras.txt --taxonomy uparse/NewTaxa.missed.txt; 
mv otu_table.biom uparse/otu_table_missed.biom;
mv reads2otus.txt uparse/reads2otus.missed.txt;
perl ~mcwong/GitRepos/mcwong/ITS_workflows/make_uc_for_unmapped.pl uparse/missed.clean.uc > uparse/missed.unmapped.uc
mv NewTaxa.unmapped.txt uparse/NewTaxa.unmapped.txt;
perl ~mcwong/GitRepos/mcwong/16S_workflows/resolveIterativeUparse.pl uparse/cluster_* uparse/singletons2otus.uc uparse/missed.unmapped.uc   --derep uparse/derep.uc --chimeras uparse/chimeras.txt --taxonomy uparse/NewTaxa.unmapped.txt;
mv reads2otus.txt uparse/reads2otus.unmapped.txt;
mv otu_table.biom uparse/otu_table_unmapped.biom;
biom summarize-table -i uparse/otu_table_97.biom -o uparse/stats.otu_table_97.txt;
biom summarize-table -i uparse/otu_table_missed.biom -o uparse/stats.otu_table_missed.txt;
biom summarize-table -i uparse/otu_table_unmapped.biom -o uparse/stats.otu_table_unmapped.txt;
cat split_libraries/seqs.fna | grep "^>" | cut -f1 -d "_" | cut -f2 -d ">" | sort | uniq -c > uparse/Stats.MergedReads.txt;
cat uparse/stats.otu_table_97.txt | tail -n +17 | sed 's/^ //g' | sed -re 's/: /\t/g' | sed 's/\.0$//g' > uparse/Stats.MappedReads.97.txt;
cat uparse/stats.otu_table_missed.txt | tail -n +17 | sed 's/^ //g' | sed -re 's/: /\t/g' | sed 's/\.0$//g' > uparse/Stats.MappedReads.missed.txt;
cat uparse/stats.otu_table_unmapped.txt | tail -n +17 | sed 's/^ //g' | sed -re 's/: /\t/g' | sed 's/\.0$//g' > uparse/Stats.MappedReads.unmapped.txt;
perl ~mcwong/GitRepos/mcwong/ITS_workflows/StatsComparisonMergedVsMapped.pl uparse/Stats.MergedReads.txt uparse/Stats.MappedReads.97.txt uparse/Stats.MappedReads.missed.txt uparse/Stats.MappedReads.unmapped.txt > uparse/Stats.Combined.txt;
cat uparse/missed.centroids.uc | grep "^N" | cut -f9 | ~mcwong/getSeq uparse/missed.fa  > uparse/OTUsFailedToMap2ndTime.fa;
