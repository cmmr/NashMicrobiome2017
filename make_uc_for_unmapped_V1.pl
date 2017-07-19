#!/usr/bin/env perl
#
use warnings;
use strict;



my $raw = $ARGV[0];
my $picked = $ARGV[1];
my $low = $ARGV[2];
my $missed = $ARGV[3];



my $reads = {};
open IN, "$raw";

my @raw;
open STATSRAW, ">tempraw";
open STATSPICK, ">temppick";
open STATSLOW, ">templow";
open STATSMISSED, ">tempmissed";
while (my $line = <IN>) {
	chomp $line;
	$line =~ s/^\s+//g;
	my @parts = split /\s+/, $line;
	$reads->{$parts[1]}->{'raw'} = $parts[0];
	print STATSRAW "$parts[0]\n";
}
close STATSRAW;
close IN;



open IN, "$picked";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
#	print "sample = $parts[0]\n";
	$reads->{$parts[0]}->{'picked'} = $parts[1];
	print STATSPICK "$parts[1]\n";
}
close STATSPICK;
close IN;
open IN, "$low";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	$reads->{$parts[0]}->{'low'} = $parts[1];
	print STATSLOW "$parts[1]\n";
}
close STATSLOW;
close IN;
open IN, "$missed";
while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	$reads->{$parts[0]}->{'missed'} = $parts[1];
	print STATSMISSED "$parts[1]\n";
}
close IN;
close STATSMISSED;
foreach my $sample (keys %{$reads}) {
	unless (exists $reads->{$sample}->{'picked'}) {
		$reads->{$sample}->{'picked'} = 0
	}
	unless (exists $reads->{$sample}->{'raw'}) {
		$reads->{$sample}->{'raw'} = 0;
	}
	unless(exists $reads->{$sample}->{'low'}) {
		$reads->{$sample}->{'low'} = 0;
	}
	unless(exists $reads->{$sample}->{'missed'}) {
		$reads->{$sample}->{'missed'} = 0;
	}

}


print "SampleID\tMapped\tLow Identity Mapped\tUnmapped\tRemainder=Chimeras+Singletons\n";
foreach my $sample (sort {$reads->{$a}->{'raw'} <=> $reads->{$b}->{'raw'}} keys %{$reads}) {
	my $picked = $reads->{$sample}->{'picked'};
	my $low = $reads->{$sample}->{'low'};
	my $unmapped = $reads->{$sample}->{'missed'};
	my $diff = $reads->{$sample}->{'raw'} - $picked - $low - $unmapped;
	print "$sample\t$picked\t$low\t$unmapped\t$diff\n";
}

print "\n\n";
print "Raw Stats:\n";
my $capture = `cat tempraw | ~mcwong/listStats.pl`;
print "$capture";
`rm tempraw`;
print "\n\n";
print "Mapped Stats:\n";
$capture = `cat temppick | ~mcwong/listStats.pl`;
print "$capture";
`rm temppick`;
print "\n\n";
print "Low Identity Hits:\n";
$capture = `cat templow | ~mcwong/listStats.pl`;
print "$capture";
`rm templow`;
print "\n\n";
print "Unmapped Reads:\n";
$capture = `cat tempmissed | ~mcwong/listStats.pl`;
print "$capture";
`rm tempmissed`;
[auchtung@cmmr-login1 ITS_workflows]$ cat make_uc_for_unmapped.pl 
#!/usr/bin/env perl
#
#
use warnings;
use strict;


open IN, "$ARGV[0]";

my $unknownCount = 0;

while (my $line = <IN>) {
	chomp $line;
	my @parts = split /\t/, $line;
	if ($parts[0] =~ m/^H/) {
		$parts[9] = "*";
	} else {
		$parts[9] = "UNMAPPED_OTU_$unknownCount";
		$parts[0] = "H";
		$unknownCount++;
	}
	my $newline = join "\t", @parts;
	print "$newline\n";

}

close IN;

open OUT, ">NewTaxa.unmapped.txt";

for (my $i = 0; $i < $unknownCount; $i++) {

	print OUT "UNMAPPED_OTU_$i\tk__Unmapped; p__Unmapped; c__Unmapped; o__Unmapped; f__Unmapped; g__Unmapped; s__Unmapped\n";
}