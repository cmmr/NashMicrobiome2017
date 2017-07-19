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
[auchtung@cmmr-login1 ITS_workflows]$ cat cleanHitsTableITS.pl 
#!/usr/bin/env perl
#
use warnings;
use strict;
#use DataBrowser qw(browse browseErr);

open IN, "$ARGV[0]";
my $taxa = {};
open TAXA, "$ARGV[1]";
my $maxslice = 0;
while (my $line = <TAXA>) {
	chomp $line;
	$line =~ s/; /;/g;
	my @parts = split /\t/, $line;
	my @taxaSlice = split /;/, $parts[1];
	unless ($maxslice) {
		$maxslice = scalar(@taxaSlice) - 1;
	}
	$taxa->{$parts[0]} = $parts[1];

}
close TAXA;

#browse($taxa);
my $otus = {};
my $badUsed = {};
#grabs all matches for an otu
while (my $line = <IN>) {

	chomp $line;
	my @parts = split /\t/, $line;
	my $queryLength = $parts[2];
	my $id = $parts[3];
	my $matchString = $parts[7];
	#this assumes the first match amount will be the biggest. Can be corrected if edge cases appear.
	$matchString =~ m/(\d+)M/;
	my $matchLength = $1;
	unless ($matchLength) {
		print "$line\n";;
		next;
	}
	my $otu = $parts[8];
	my $taxaId = $parts[9];
	#compensating for quirks in the .uc file.
	if (($id eq "*" or $taxaId eq "*") and not exists $badUsed->{$otu}) {
		print "$line\n";
		$badUsed->{$otu} = 1;
		next;
	}
	if (($id eq "*" or $taxaId eq "*") and exists $badUsed->{$otu}) {
		next;
	}
	$otus->{$otu}->{$id}->{$matchLength}->{$taxaId} = $taxa->{$taxaId};
}
close IN;
my $candidates = {};
my $topId = {};
my $topIdName = {};

#cleans up the otu hits into the candidate top hits only.
#-top-hits-only in usearch has weird quirks when it comes to leading or terminal deletions in the reference. IE. I3290M might be as good as 293M
#we maximize both identity and length of identity
#browse($otus);
foreach my $otu (keys %{$otus}) {
	my $maxid = 0;
	my $matchAmount = 0;
	foreach my $id (keys %{$otus->{$otu}}) {
		my $maxMatch = 0;
		if ($id < $maxid) {
			delete $otus->{$otu}->{$id};
			next;
		} 
		if ($id > $maxid) {	
			unless ($maxid == 0) {
				delete $otus->{$otu}->{$maxid};
			}
			$maxid = $id;
		}
		foreach my $matchLength (keys %{$otus->{$otu}->{$id}}) {
			if ($matchLength < $maxMatch) {
				delete $otus->{$otu}->{$id}->{$matchLength};
				next;
			}
			if ($matchLength > $maxMatch) {
				unless ($maxMatch == 0) {
					delete $otus->{$otu}->{$id}->{$maxMatch};
				}
				$maxMatch = $matchLength;
			}
		}
		$matchAmount = $maxMatch;
	}
	$candidates->{$otu} = $otus->{$otu}->{$maxid}->{$matchAmount};
	$topId->{$otu} = $maxid;
}
#browse($candidates);
my $taxaCount = 0;
open TAXA, ">NewTaxa.txt";
#print "before\n";
foreach my $otu (keys %{$candidates}) {
	#print "here\n";
	#print "START $otu\n";
	if (exists $badUsed->{$otu}) {
		next;
	}
	$taxaCount++;
	#browse($candidates->{$otu});
	my @candidateOtus = keys %{$candidates->{$otu}};
	if (scalar(@candidateOtus) == 1) {
		print "H\t*\t*\t$topId->{$otu}\t+\t0\t0\t*\t$otu\t$candidateOtus[0]\t1\n";	
		$candidates->{$otu}->{$candidateOtus[0]} =~ s/;/; /g;
		print TAXA "$candidateOtus[0]\t$candidates->{$otu}->{$candidateOtus[0]}\n";
		next;
	}
	my $levels = [];
	my $names = {};
	#separate each candidate otu into identifiable parts, then assign those parts into a candidate new name hash
	#simultaneously count each occurence of a clade per level
	#names are used only if the current clade is identifable. Subsequent levels are always stripped off.
	foreach my $taxaId (@candidateOtus) {
#		print "taxaId = $taxaId\n";
#		print "candidates = $candidates->{$otu}->{$taxaId}\n";
		my @parts = split /;/, $candidates->{$otu}->{$taxaId};
		for (my $i = 0; $i <= $maxslice; $i++) {
			if ($parts[$i] =~ m/unidentified/ or $parts[$i] =~ m/uncultured/ or $parts[$i] =~ m/ sp\.\s?/) {
				next;
			}
			my @slice = @parts[0..$i];
			my $slice = join ";", @slice;
#			print "slice = $slice\n";
			$levels->[$i]->{$slice}++;
			$names->{$slice}->{$taxaId} = 1;
		}
	}
#	browse($levels);
	my $found = 0;
	my $taxonomy = "";
	my $highestLevel = 0;
	my $lastName= "";
	#starting at the highest classification of taxonomy, determine if there is at least 80% agreement against the total candidates in that clade. (This could change to total candidates).
	#if so, use that level of classification.
	#the key used to find the best matching accession is the trimmed name. Since $names stores each truncated name's accession, we have a pool to pick from
	for (my $i = $maxslice; $i >= 0; $i--) {
		my $total = 0;
		if ($found) {
			last;
		}
#		print "i = $i\n";
		foreach my $candidate (keys %{$levels->[$i]}) {
			$total += $levels->[$i]->{$candidate};
		}
#		print "total = $total\n";
		foreach my $candidate (keys %{$levels->[$i]}) {
			if (($levels->[$i]->{$candidate} / $total) > .8) {
				my @parts = split ";", $candidate;
				$lastName = $parts[$i];
				$lastName =~ s/.*__//g;
				$found = 1;
				$taxonomy = $candidate;
				$highestLevel = $i;
			}
		}
	}
	my @possibleTaxaNames = sort keys %{$names->{$taxonomy}};
	my $taxaName = "";
	if (scalar(@possibleTaxaNames)) {
		$taxaName = $possibleTaxaNames[0];
	} else {
		$taxaName = $candidateOtus[0];
	}
	unless($taxonomy) {
		print "H\t*\t*\t*\t+\t0\t0\t*\t$otu\t*\n";
		print TAXA "$taxaName\tk__ambiguous; p__ambiguous; c__ambiguous; o__ambiguous; f__ambiguous; g__ambiguous; s__ ambiguous\n";
		next;
	}
#	print "highest level = $highestLevel\n";
#	browse($names);
#	print "otu = $otu\n";
#	browse($levels);
	#this part fills in the highest level of confidence allowed by the picking method, and appends <clade> sp. 
	my @taxaLevel = qw(k p c o f g s);
	for (my $i = ($highestLevel + 1); $i <= $maxslice; $i++) {
		$taxonomy = $taxonomy . ";$taxaLevel[$i]__$lastName sp.";
	}
	$taxonomy =~ s/;/; /g;
	#my @newLine;
	print "H\t*\t*\t$topId->{$otu}\t+\t0\t0\t*\t$otu\t$taxaName\n"; 
	print TAXA "$taxaName\t$taxonomy\n";
	#print "$taxaName\t$taxonomy\n";
	#print  "END $otu\n\n";
}
close TAXA;	



#browse($candidates);