#!/usr/bin/env perl
#
use warnings;
use strict;
use DataBrowser qw(browse browseErr);

my $read1 = $ARGV[0];
my $read2 = $ARGV[1];
my $barcodeList = $ARGV[2];

open IN, "$barcodeList";


my $barcodes = {};
my $rcbarcodes = {};
my $bclengths = {};
my $ambig = {};
$ambig->{"R"} = [qw(A G)];
$ambig->{"Y"} = [qw(C T)];
$ambig->{"M"} = [qw(A C)];
$ambig->{"K"} = [qw(G T)];
$ambig->{"S"} = [qw(C G)];
$ambig->{"W"} = [qw(A T)];
$ambig->{"H"} = [qw(A C T)];
$ambig->{"B"} = [qw(C G T)];
$ambig->{"V"} = [qw(A C G)];
$ambig->{"D"} = [qw(A G T)];
$ambig->{"N"} = [qw(A C G T)];


sub makeTree {
	my $point = $_[0];
	my $primer = $_[1];
	my $name = $_[2];
	my @front = split //, $primer;
        my $prevPoint;
        my $lastLetter;
        foreach my $letter (@front) {
                my $hash = {};
                if (exists $ambig->{$letter}) {
                        foreach my $amLet (@{$ambig->{$letter}}) {
                                $point->{$amLet} = $hash;
                        }
                        $point = $hash;
                        next;
                } elsif (not (exists $point->{$letter})) {
                        $point->{$letter} = $hash;
                }
                $prevPoint = $point;
                $lastLetter = $letter;
                $point = $point->{$letter};
        }
	$prevPoint->{$lastLetter} = $name;

}
my %fh;
my $used = {};
while (my $line = <IN>) {
	chomp $line;
	my @fields = split /\t/, $line;
	unless (scalar(@fields) == 3) {
		warn "line: $line is not properly formatted\n";
		next;
	}
	my $name = $fields[0];
	if (exists ($used->{$name})) {
		die "repeated name\n";
	}
	$used->{$name} = 1;
	open ($fh{"$name.1"}, ">$name.1.fq");
	open ($fh{"$name.2"}, ">$name.2.fq");
	makeTree($barcodes, $fields[1], $fields[1]);
	$bclengths->{$name}->{"1"} = length($fields[1]);
	$bclengths->{$name}->{"2"} = length($fields[2]);
	unless(exists $rcbarcodes->{$fields[1]}) {
		$rcbarcodes->{$fields[1]} = {};
	}
	makeTree($rcbarcodes->{$fields[1]}, $fields[2], $name);
}
#browse($bclengths);
#browse($barcodes);
#browse($rcbarcodes);
#die;
open R1, "zcat $read1 |";
open R2, "zcat $read2 |";

my $count = 0;

sub resolveTree {
	my $seq = reverse($_[1]);
        my $point = $_[0];
	my $start = $_[2];
	my $mm = 0;
	if ($start) {
		$point = $point->{$start};
	}
        my $end = "";
        while(1) {
                my $letter = chop $seq;
		#print "letter = $letter\n";
                if (exists $point->{$letter} and ref($point->{$letter}) eq "HASH" ) {
                        $point = $point->{$letter};
                        next;
                } elsif (exists $point->{$letter}) {
                        $end = $point->{$letter};
                        last;
                } elsif ($mm > 2) {
		#	print "end\n";
                        last;
                } else {
			last;
			$mm += 1;
			my @keys = keys %{$point};
			my $key = $keys[0];
	#		unless ($point->{$key}) {
	#			browseErr($point)
	#		}	
	#		print STDERR "top key = $key\n";
			if (($mm <= 2) and ref($point->{$key}) eq "HASH") {
		#		print "here as hash\n";
				$point = $point->{$key};
			} elsif ($mm <= 2) {
		#		print "not as hash\n";
				$end = $point->{$key};
				last	
			}
		}
        }
	return($end);

}
my $multi = 0;
my $nobc = 0;
my $sampleCounts = {};
while (my $head1 = <R1>) {
	my $seq1 = <R1>;
	my $space1 = <R1>;
	my $qual1 = <R1>;
	my $head2 = <R2>;
	my $seq2 = <R2>;
	my $space2 = <R2>;
	my $qual2 = <R2>;
	#$count++;
	#if ($count % 10000 == 0) {
	#	print STDERR "processed $count reads\n";
	#}

	my $bc1 = resolveTree($barcodes, $seq1);
	#print "seq1= $seq1\n";
	#print "bc1 = $bc1\n";
	unless ($bc1) {
		next;
	}
	#print "bc1 = $bc1\n";
	my $sample = resolveTree($rcbarcodes, $seq2, $bc1);
	#print "bc2 = $sample\n";
	unless ($sample) {
		next;
	}
	#print "bc1 = $bc1\n";
	#print "sample = $sample\n";
	my $length1 = length($seq1);
	my $length2 = length($seq2);
	#print "origseq1 = $seq1";
	my $cut1 = $bclengths->{$sample}->{"1"};
	my $cut2 = $bclengths->{$sample}->{"2"};
	$seq1 = substr($seq1, $cut1, $length1 - $cut1);
	#print "new seq1 = $seq1";
	#print "origseq2 = $seq2";
	$seq2 = substr($seq2, $cut2, $length2 - $cut2);
	$qual1 = substr($qual1, $cut1, $length1 - $cut1);
	$qual2 = substr($qual2, $cut2, $length2 - $cut2);
	#print "new seq2 = $seq2";
	#print "\n";
	print {$fh{"$sample.1"}} "$head1$seq1$space1$qual1";
	print {$fh{"$sample.2"}} "$head2$seq2$space2$qual2";
	#print TEST1 "$head1$seq1$space1$qual1";
	#print TEST2 "$head2$seq2$space2$qual2";
	$sampleCounts->{$sample}++;
}
browseErr($sampleCounts);
#print STDERR "Total raw reads: $count\n";
#print STDERR "Total no barcode: $nobc\n";
#print STDERR "Total multi barcode: $multi\n";