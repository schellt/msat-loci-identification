#!/usr/bin/env perl
use warnings;
# Author: Thomas Thiel
# Adapted: Tilman Schell
# Program name: primer3_in.pl
# Description: creates a PRIMER3 input file based on SSR search results

open (IN,"<$ARGV[0]") || die ("\nError: Couldn't open misa.pl results file (*.misa) !\n\n");

my $filename = $ARGV[0];
$filename =~ s/\.misa_filter//;
open (SRC,"<$filename") || die ("\nError: Couldn't open source file containing original FASTA sequences !\n\n");
open (OUT,">$filename.p3in");

undef $/;
$in = <IN>;
study $in;

$/= ">";

my $count = 0;
while (<SRC>)
  {
  next unless (my ($id,$seq) = /(.*?)\n(.*)/s);
  $seq =~ s/[\d\s>]//g;#remove digits, spaces, line breaks,...
  while ($in =~ /$id\t(\d+)\t\S+\t\S+\t(\d+)\t(\d+)/g)
    {
    my ($ssr_nr,$size,$start) = ($1,$2,$3);
    $count++;
    #print OUT "PRIMER_SEQUENCE_ID=$id"."_$ssr_nr\nSEQUENCE=$seq\n";
    print OUT "PRIMER_SEQUENCE_ID=$id"."_$ssr_nr\nSEQUENCE_TEMPLATE=$seq\n";
    print OUT "PRIMER_PRODUCT_SIZE_RANGE=100-300\n";
    print OUT "PRIMER_PRODUCT_OPT_SIZE=150\n";
#    print OUT "PRIMER_NUM_RETURN=1\n";
    print OUT "PRIMER_MAX_END_STABILITY=9.0\n";
    print OUT "PRIMER_MAX_LIBRARY_MISPRIMING=12.00\n";
    print OUT "PRIMER_PAIR_MAX_LIBRARY_MISPRIMING=24.00\n";
    print OUT "PRIMER_MIN_SIZE=18\n";
    print OUT "PRIMER_OPT_SIZE=21\n";
    print OUT "PRIMER_MAX_SIZE=23\n";
    print OUT "PRIMER_MIN_TM=50\n";
    print OUT "PRIMER_OPT_TM=55\n";
    print OUT "PRIMER_MAX_TM=60\n";
    print OUT "PRIMER_PAIR_MAX_DIFF_TM=8.0\n";
    print OUT "PRIMER_MIN_GC=30.0\n";
    print OUT "PRIMER_OPT_GC_PERCENT=50.0\n";
    print OUT "PRIMER_MAX_GC=70.0\n";
    print OUT "PRIMER_MAX_SELF_ANY=8.00\n";
    print OUT "PRIMER_MAX_SELF_END=3.00\n";
    print OUT "PRIMER_MAX_NS_ACCEPTED=0\n";
    print OUT "PRIMER_MAX_POLY_X=5\n";
    print OUT "PRIMER_OUTSIDE_PENALTY=0\n";
    print OUT "PRIMER_GC_CLAMP=0\n";
    print OUT "PRIMER_SALT_MONOVALENT=50.0\n";
    print OUT "PRIMER_DNA_CONC=50.0\n";
    print OUT "PRIMER_LIBERAL_BASE=1\n";
    print OUT "SEQUENCE_TARGET=",$start-3,",",$size+6,"\n";
    print OUT "=\n"
    };
  };
print "\n$count records created.\n";
