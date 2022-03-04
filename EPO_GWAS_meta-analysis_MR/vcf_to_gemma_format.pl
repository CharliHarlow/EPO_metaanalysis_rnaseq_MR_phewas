#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Long;

## Converting from vcf format to gemma format 

my $_vcf  = undef;
my $_out  = undef;
my $_line = undef;

GetOptions
(
  "vcf=s" => \$_vcf,
  "out=s" => \$_out,
);

open VCF, "gunzip -c $_vcf |" or die "Unable to open vcf for reading\n";
open OUT, ">" . $_out or die "Unable to open output file\n";

while ($_line = <VCF>) {
    if ($_line =~ /#CHROM/) {
        last;
    }
}

while ($_line = <VCF>) {
    chomp $_line;
    my @_data = split(/\t/, $_line);
    my $_stringout =  $_data[2] . "," . $_data[3] . "," . $_data[4];
    for (my $_i = 9; $_i < @_data; $_i++) {
        my @_geno_data = split(":", $_data[$_i]);
        $_stringout .= "," . $_geno_data[1];
    }
    $_stringout .= "\n";
    print OUT $_stringout;
}

close OUT;
close VCF;
