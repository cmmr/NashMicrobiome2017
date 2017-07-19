#!/usr/bin/env perl
#
use warnings;
use strict;
#use Text::Levenshtein::XS qw(distance);
use Text::LevenshteinXS qw(distance);
use DataBrowser qw(browse);
open R1, "$ARGV[0]";
open R2, "$ARGV[1]";
open MER, "$ARGV[2]";



my $read1 = {};
my $read2 = {};

while (my $head1 = <R1>) {
	my $head2 = <R2>;
	my $seq1 = <R1>;
	my $seq2 = <R2>;
	my $space1 = <R1>;
	my $space2 = <R2>;
	my $qual1 = <R1>;
	my $qual2 = <R2>;
	chomp $head1;
	chomp $head2;
	chomp $seq1;
	chomp $seq2;
	chomp $space1;
	chomp $space2;
	chomp $qual1;
	chomp $qual2;
	$head1 =~ s/ .*//g;
	$head2 =~ s/ .*//g;
	my $prefix1 = substr($seq1, 0, 20);
	my $prefix2 = substr($seq2, 0, 20);
	$read1->{$head1} = $prefix1;
	$read2->{$head2} = $prefix2;
}
close R1;
close R2;
my $trimmedFront = 0;
my $trimmedBack = 0;
my $totalReads = 0;
while (my $head = <MER>) {
	$totalReads++;
	chomp $head;
	my $seq = <MER>;
	chomp $seq;
	my $space = <MER>;
	chomp $space;
	my $qual = <MER>;
	chomp $qual;
	my $key = $head;
	$key =~ s/ .*//g;
	my $front = substr($seq, 0, 20);
	my $back = substr($seq, length($seq) - 20, 20);
	$back =~ tr/ACGT/TGCA/;
	$back = reverse($back);
	my $origFront = $read1->{$key};
	my $origBack = $read2->{$key};
	if (($front ne $origFront) or ($back ne $origBack)) {
		#print "$head\n";
		#print "orig = \n$seq\n$qual\n";
		#print "origFront = $origFront\n";
		#print "nowFront  = $front\n";
		#print "origBack = $origBack\n";
		#print "nowBack = $back\n";
		my $offsetF = -1;
		my $offsetB = -1;
		my $lowestDistanceF = -1;
		my $lowestDistanceB = -1;
		for (my $i = 0; $i < length($seq) - 20; $i++) {
			my $newFront = substr($seq, $i, 20);
			if ($newFront eq $origFront) {
				$offsetF = $i;
				last
			}
			my $distance = distance($newFront,$origFront);
			if ($distance < 3 and $lowestDistanceF == -1) {
				$lowestDistanceF = $distance;
				$offsetF = $i;
				next;
			} elsif ($distance < 3) {
				if ($distance > $lowestDistanceF) {
					last;
				} else {
					$lowestDistanceF = $distance;
					$offsetF = $i;
					next;
				}
			}
		}
		if ($offsetF == -1) {
			#print "read was not recoverable for read 1\n";
		} elsif ($offsetF == 0) {
			#print "read was correct on the front\n";
		} else {
			$trimmedFront++;
			$seq = substr($seq, $offsetF, (length($seq) - $offsetF));
			$qual = substr($qual, $offsetF, (length($qual) - $offsetF));
			#print "new1 =\n$seq\n$qual\n";
			#print "\n";
		}
		$seq = reverse($seq);
		$seq =~ tr/ACGT/TGCA/;
		$qual = reverse($qual);
		for (my $i = 0; $i < length($seq) - 20; $i++) {
			my $newBack = substr($seq, $i, 20);
			if ($newBack eq $origBack) {
				$offsetB = $i;
				last;
			}
			my $distance = distance($newBack,$origBack);
			if ($distance < 3 and $lowestDistanceB == -1) {
				$lowestDistanceB = $distance;
				$offsetB = $i;
				next;
			} elsif ($distance < 3) {
				if ($distance > $lowestDistanceB) {
					last;
				} else {
					$lowestDistanceB = $distance;
					$offsetB = $i;
					next;
				}
			}
		}
		if ($offsetB == -1) {
			#print "read was not recoverable for read 2\n";
		} elsif ($offsetB == 0) {
			#print "read was correct on the back\n";
		} else {
			$trimmedBack++;
			$seq = substr($seq, $offsetB, (length($seq) - $offsetB));
			$qual = substr($qual, $offsetB, (length($qual) - $offsetB));
			#print "new2 =\n$seq\n$qual\n";
		}
		$seq =~ tr/ACGT/TGCA/;
		$seq = reverse($seq);
		#print "$key\n";
		#print "front  = $front\n";
		#print "back   = $back\n";
		#print "\n";
		#print "origF  = $origFront\n";
		#print "origB  = $origBack\n";
		#print "\n\n";
	}
	print "$head\n";
	print "$seq\n";
	print "$space\n";
	print "$qual\n";
}
close MER;
print STDERR "trimmed $trimmedFront/$totalReads reads from the front\n";
print STDERR "trimmed $trimmedBack/$totalReads reads from the back\n";
print STDERR "\n";