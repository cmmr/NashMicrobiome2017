#!/usr/bin/env perl
#
my $in = $ARGV[0];

my $barcodes = $ARGV[1];

unless($barcodes && $in) {
        die "need both source and barcode file\nusage:barcoderecover.pl <source> <barcodes>\n";
}


my $cmd = "cat $barcodes";
if ($barcodes =~ m/\.gz$/) {
        $cmd = "zcat $barcodes";

} elsif ($barcodes =~ m/\.bz2$/) {
        $cmd = "pbzip2 -c -d -p16 $barcodes";
}

open BAR, "$cmd |";

my %barcodes;
my $count = 0;
while (my $header = <BAR>) {
        my $seq = <BAR>;
        my $space = <BAR>;
        my $qual = <BAR>;

        my @fields = split /\s/, $header;
        
        my $string = "$header$seq$space$qual";
        $barcodes{$fields[0]} = $string;
#       if ($count < 10) {
#               print STDERR "added $fields[0]\n";
#               $count++;
#       }

}

close BAR;
$cmd = "cat $in";

if ($in =~ m/\.gz$/) {
        $cmd = "zcat $in";
} elsif ($in =~ m/\.bz2$/) {
        $cmd = "pbzip2 -c -d -p8 $in";
}

open MER, "$cmd |";
#print STDERR "Checking source file now";

while (my $header = <MER>) {
        my $seq = <MER>;
        my $space = <MER>;
        my $qual = <MER>;

        my @fields = split /\s/, $header;
        if (exists $barcodes{$fields[0]}) {
                print "$barcodes{$fields[0]}";
        } else {
                print STDERR "$fields[0] is in source file, but not in barcodes file\n";
                die "$fields[0] is not found in barcodes file\n";

        }


}