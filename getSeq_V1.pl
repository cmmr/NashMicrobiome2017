#!/usr/bin/env perl


use warnings;
use strict;
use FAlite;

my $reverse = 0;
if ($ARGV[0] eq "-v") {
	$reverse = 1;
	shift @ARGV;
}
my $in = $ARGV[0];
my $seqhead = $ARGV[1];
my $fromSTDIN = 0;
my %headers;
unless ($seqhead) {
	open HEAD, "-";
	while (<HEAD>) {
		my $line = $_;
		my @fields = split(/\s+/, $line);
		my $accession = $fields[0];
		chomp $line;
		$headers{$line} = 1;
		$headers{$accession} = 1;
	}
} else {
	$headers{$seqhead} = 1;
}


open IN, "$in";
my %seqs;

my $fasta_file = new FAlite(\*IN); # or any other filehandle
while(my $entry = $fasta_file->nextEntry) {
	my $head = $entry->def;
	$head =~ s/>//g;
	my $fullName = "";
	#print "head = $head\n";
	if ($head =~ m/(\S+);size=/) {
		$fullName = $1;
	}
	#print "fullName = $fullName\n";
	#my @fields = split(/\s+/, $head);
	#my $fullName = $fields[0];
	unless ($reverse) {
		if (exists ($headers{$head}) or exists ($headers{$fullName})) {
			#if (exists ($headers{$fullName})) {
			my $seq = uc($entry->seq);
			print ">$head\n";
			print "$seq\n";
		}
	} else {
		unless (exists ($headers{$head}) or exists ($headers{$fullName})) {
			my $seq = uc($entry->seq);
			print ">$head\n";
			print "$seq\n";
		}	
	}
}
