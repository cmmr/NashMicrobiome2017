#!/usr/bin/env perl
#
#
use warnings;
use strict;
use Getopt::Long;
use DataBrowser qw(browse);
use JSON::XS;
use JSON;
use Time::HiRes qw(gettimeofday);

my $derep;
my $chimeras;
my $uchime;
my $taxonomy;
my $all;
my $writePerc;
my $otuFile;
my $pipelineName;
GetOptions ("perc" => \$writePerc, "pipeline=s" => \$pipelineName, "derep=s" => \$derep, "chimeras=s" => \$chimeras, "uchime=s" => \$uchime, "taxonomy=s" => \$taxonomy, "otus=s" => \$otuFile, "all" => \$all) or die "Can't GetOpts\n";


unless ($derep) {
	die "needs --derep derep.uc --chimeras chimeras.txt in command line\n";
}
unless ($taxonomy) {
	die "needs --taxonomy in the command line\n";
}
if ($writePerc and not $otuFile) {
	die "--perc needs --otus\n";
}
my $pipelinePath = $ENV{"AMPLICON_PIPELINES"};
my $pipelineFile = "";
if ($pipelinePath and $pipelineName and -d $pipelinePath and -e "$pipelinePath/$pipelineName") {
	$pipelineFile = "$pipelinePath/$pipelineName";
}
#print "pipelinePath = $pipelineFile";
#die;
my %tree;
my %alone;
foreach my $file (@ARGV) {
	if ($file =~ m/cluster_(\d\.\d)\./) {
		my $radius = $1;
		print STDERR "Reading clustering $file\n";
		open IN, "$file" or die "can't open $file\n";
		my $minsize = "inf";
		while (my $line = <IN>) {
			chomp $line;
			my @fields = split /\t/, $line;
			unless (scalar(@fields) == 10) {
				last;
			}
			if ($line =~ m/^C/) {
				last
			}
			unless ($fields[8]) {
				next;
			}
			my $clusterSize;
			if ($fields[8] =~ m/size=(\d+);/) {
				$clusterSize = $1;
			} else {
				next;
			}
			if ($clusterSize < $minsize) {
				$minsize = $clusterSize
			} elsif ($clusterSize > $minsize) {
				next
			}
			$fields[8] =~ s/;.*//g;
			$fields[9] =~ s/;.*//g;
			if ($line =~ m/^S/) {
				$alone{$fields[8]} = "*";
			}
			unless ($line =~ m/^H/) {
				next;
			}
			my $limit = 100 - $radius;
			my $overlap = $fields[3];
			if ($fields[3] < $limit) {
				next;
			}
			if (exists $alone{$fields[8]}) {
				delete ($alone{$fields[8]});
			}
			if ($fields[8] eq $fields[9]) {
				die "$line broke\n";
			}
			if (not exists $tree{$fields[9]} and not exists $alone{$fields[9]}) {
				print STDERR "Dead node $fields[8]\n";
				$tree{$fields[8]} = "dead";
				next;
			}
			$tree{$fields[8]} = $fields[9];
		}
		close IN;
	} else {
		print STDERR "Reading standard $file\n";
		open IN, "$file" or die "can't open $file\n";
		while (my $line = <IN>) {
			chomp $line;
			unless ($line =~ m/^H/) {
				next
			}
			my @fields = split /\t/, $line;
#			if ($file =~ m/clean/) {
#				print "line = $line\n";
#				print "fields 8 = $fields[8]\n";
#				print "fields 9 = $fields[9]\n";
#			}
			$fields[8] =~ s/;.*//g;
			$fields[9] =~ s/;.*//g;
			if ($fields[8] eq $fields[9]) {
				die "$line broke\n";
			}
			$tree{$fields[8]} = $fields[9];
		}
		close IN;
	}

}
close IN;
if ($chimeras) {
	print STDERR "Reading chimeras $chimeras\n";
	open IN, "$chimeras";
	while (my $line = <IN>) {
		chomp $line;
		my @fields = split /\t/, $line;
		$fields[0] =~ s/;.*//g;
		$tree{$fields[0]} = "denovo chimera";
	}
	close IN;

}
if ($uchime) {
	print STDERR "Reading uchime $uchime\n";
	open IN, "$uchime";
	while (my $line = <IN>) {
		chomp $line;
		my @fields = split /\t/, $line;
		if ($fields[-1] eq "Y") {
			$fields[1] =~ s/;.*//g;
			$tree{$fields[1]} = "uchime chimera";
		}
	}
	close IN;
}
open IN, "$derep";
print STDERR "Reading derep $derep\n";
while (my $line = <IN>) {
	chomp $line;
	my @fields = split /\t/, $line;
	if ($line =~ m/^S/) {
		unless (exists $tree{$fields[8]}) {
			$tree{$fields[8]} = $fields[9];
		}
		next
	}
	unless ($line =~ m/^H/) {
		next
	}
	unless (exists $tree{$fields[9]}) {
		$tree{$fields[8]} = "*";
		next
	}
	$tree{$fields[8]} = $fields[9];

}
close IN;

print STDERR "Resolving Tree\n";

foreach my $alone (keys %alone) {
	if (exists $tree{$alone}) {
		next
	}
	$tree{$alone} = $alone{$alone};

}
#browse(\%tree);


print STDERR "Tracing Lineage\n";
open OUT, ">readLineageTrace.out";
foreach my $leaf (keys %tree) {
	my @trace;
	my $origLeaf = $leaf;
	push @trace, $leaf;
	push @trace, $tree{$leaf};
	#print "Original Parent = $tree{$leaf}\n";
	my $traceCount = 0;
	while (1) {
		$traceCount++;
		if ($traceCount == 100) {
			print STDERR "Over 100 traces on $origLeaf, exiting trace\n";
			last;
		}
		if (exists ($tree{$tree{$leaf}})) {
			my $newParent = $tree{$tree{$leaf}};
			#print "NewParent = $newParent\n";
			push @trace, $newParent;
			$tree{$leaf} = $newParent;
		} else {
			my $trace = join " => ", @trace;
			print OUT "$trace\n";
			last
		}
	}
}
close OUT;
print STDERR "Printing reads2otus.txt\n";
my $otus = {};
my $samples = {};
open OUT, ">reads2otus.txt";
my $readQR = qr/_\d+.*$/;
foreach my $node (keys %tree) {
	print OUT "$node\t$tree{$node}\n";
	if ($tree{$node} eq "*" or $tree{$node} eq "denovo chimera" or $tree{$node} eq "uchime chimera" or $tree{$node} eq "dead") {
		next
	}
	my $sampleName = $node;
	$sampleName =~ s/$readQR//;
	$samples->{$sampleName} = 1;
	$otus->{$tree{$node}}->{$sampleName}++;
}
close OUT;
my $taxa = {};
print STDERR "Reading Taxa File\n";
open TAX, "$taxonomy";
while (my $line = <TAX>) {
	chomp $line;
	my @parts = split /\t/, $line;
	my @taxParts = split /; /, $parts[1];
	$taxa->{$parts[0]} = \@taxParts;
}
close TAX;

my $otuPercentages = {};
if ($writePerc and $otuFile) {
	open IN, "$otuFile";
	while (my $line = <IN>) {
		chomp $line;
		my @parts = split /\t/, $line;
		$otuPercentages->{$parts[9]} = $parts[3];
	}
}
close IN;
print STDERR "Making OTU Table\n";
my $otuTable = {};
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
$year += 1900;
$mon += 1;
$mon = sprintf("%02d", $mon);
$mday = sprintf("%02d", $mday);
my @sec = gettimeofday();
my $pipelineInfo = "";
my $pipelineJSON = "";
if ($pipelineFile) {
	open IN, "$pipelineFile";
	while (my $line = <IN>) {
		chomp $line;
		$pipelineInfo .= $line;
	}
	close IN;
	$pipelineJSON = decode_json($pipelineInfo);
} else {
	print "Warning: No pipeline or invalid pipeline defined\n";

}

$otuTable->{'date'} = "$year-$mon-$mday"."T$hour:$min:$sec.$sec[1]";
$otuTable->{'format'} = "CMMR .biom clone";
$otuTable->{'format_url'} = "NA";
if($pipelineInfo) {
	$otuTable->{'generated_by'} = $pipelineJSON;
} else {
	$otuTable->{'generated_by'} = {
		"Facility" => "CMMR"
	}	
}
$otuTable->{'id'} = "None";
$otuTable->{'matrix_element_type'} = "int";
$otuTable->{'matrix_type'} = "sparse";
$otuTable->{'type'} = "OTU table";
$otuTable->{'rows'} = ();
$otuTable->{'columns'} = ();
$otuTable->{'data'} = ();
my @otus = keys %{$otus};
my $otuCount = 0;
my @samples = keys(%{$samples});
my $sampleCount = scalar(@samples);
my $otuIdx = 0;
foreach my $otu (@otus) {
	unless (exists $taxa->{$otu} and $taxa->{$otu}) {
		die "$otu taxanomy not found\n";
	}
	#print STDERR "otuName = $otu\n";
	my $otuName = join " ", @{$taxa->{$otu}};
	#print STDERR "taxa = $otuName\n";
	if ($otuName =~ m/Francisella/  and not $all) {
		next;
	}
	my $hash = {};
	$hash->{'id'} = $otu;
	unless (exists ($taxa->{$otu})) {
		next;
	}
	$otuCount++;
	$hash->{'metadata'}->{'taxonomy'} = $taxa->{$otu};
	if ($writePerc) {
		$hash->{'metadata'}->{'taxonomy'}->[-1] .= " ($otuPercentages->{$otu}%)";
	} 
	push @{$otuTable->{'rows'}}, $hash;
	my $sampleIdx = 0;
	foreach my $sample (@samples) {
		if (exists ($otus->{$otu}->{$sample})) {
			my @array = ($otuIdx, $sampleIdx, $otus->{$otu}->{$sample});
			push @{$otuTable->{'data'}}, \@array;
		}
		$sampleIdx++;
	}
	$otuIdx++;
}
my @shape = ($otuCount, $sampleCount);
$otuTable->{'shape'} = \@shape;
foreach my $sample (@samples) {
	my $hash = {};
	$hash->{'id'} = $sample;
	$hash->{'metadata'} = "";
	undef $hash->{'metadata'};
	push @{$otuTable->{'columns'}}, $hash;
}

print STDERR "Printing OTU table\n";
open OUT, ">otu_table.biom";
my $otuText = JSON::XS->new->utf8->canonical->encode($otuTable);
print OUT "$otuText";
close OUT;