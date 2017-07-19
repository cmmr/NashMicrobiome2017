#!/usr/bin/env perl;
#
#
use warnings;
use strict;

open IN, "$ARGV[0]";
my $count = 0;
my $passed = 0;
while (my $line = <IN>) {
	$count++;
	chomp $line;
	my $seq = <IN>;
	my $space = <IN>;
	my $qual = <IN>;
	chomp $seq;
	chomp $space;
	chomp $qual;
	$line =~ m/(.*);ee=(.*);/;
	my $head = $1;
	my $ee = $2;
	my $length = length($seq);
	#if ($ee =~ m/e/) {
	#	print "here\n";
	#	print "$length\n";
	#	print "$ee\n";
	#}
	my $percErr = $ee / $length;
	if ($percErr <= .005) {
		$passed++;
		print "$head\n";
		print "$seq\n";
		print "$space\n";
		print "$qual\n";
	}
	#print "$percErr\n";

}
if ($count >= 1 ) {
	my $percPassed = sprintf "%.2f", (($passed/$count) * 100);
	print STDERR "$passed/$count passed ($percPassed%)\n";
} else {
	print STDERR "None passed\n";
}