#!/usr/bin/env perl

use strict;
use warnings;

#convert primer3 output into tsv
#usage:
#my_p3_out.pl <misa_output> <primer3_output>

my $id = "";
my $pair_id = "";
my $for_seq = "";
my $for_tm = "";
my $rev_seq = "";
my $rev_tm = "";
my $frag_size = "";

my %misa;

open (MISA, $ARGV[0]) or die "Could not open Inputfile " . $ARGV[0] . "\n";

while (my $line = <MISA>){
	chomp $line;
	my @m = split(/\t/,$line);
	$misa{$m[0] . "_" . $m[1]} = $m[3] . "\t" . $m[4];
}

close MISA;

my %results;

open (P3O, $ARGV[1]) or die "Could not open Inputfile " . $ARGV[1] . "\n";

while (my $line = <P3O>){
	chomp $line;
	if($line =~ m/^PRIMER_SEQUENCE_ID=/){
		$id = (split /=/,$line)[1];
	}
	if($line =~ m/^PRIMER_LEFT_\d_SEQUENCE=/){
		$pair_id = (split /_/,$line)[2];
		$for_seq = (split /=/,$line)[1];
		$results{$id . "_" . $pair_id} = $id . "_" . $pair_id . "\t" . $for_seq;
	}
	if($line =~ m/^PRIMER_LEFT_\d_TM=/){
		$for_tm = (split /=/,$line)[1];
		$results{$id . "_" . $pair_id} = $results{$id . "_" . $pair_id} . "\t" . $for_tm;
	}
	if($line =~ m/^PRIMER_RIGHT_\d_SEQUENCE=/){
		$rev_seq = (split /=/,$line)[1];
		$results{$id . "_" . $pair_id} = $results{$id . "_" . $pair_id} . "\t" . $rev_seq;
	}
	if($line =~ m/^PRIMER_RIGHT_\d_TM=/){
		$rev_tm = (split /=/,$line)[1];
		$results{$id . "_" . $pair_id} = $results{$id . "_" . $pair_id} . "\t" . $rev_tm;
	}
	if($line =~ m/^PRIMER_PAIR_\d_PRODUCT_SIZE=/){
		$frag_size = (split /=/,$line)[1];
		$results{$id . "_" . $pair_id} = $results{$id . "_" . $pair_id} . "\t" . $frag_size;
	}
}

close P3O;

my @sorted = sort keys(%results);

print "ID\tFOR_seq\tFOR_Tm\tREV_seq\tREV_Tm\tfragment_size\tmotif\tmotif_length\n";

for(my $i = 0; $i < scalar(@sorted); $i++){
	my $id = $sorted[$i];
	$id =~ s/_\d$//;
	my @r = split(/\t/,$results{$sorted[$i]});
	print $r[0] . "\t" . $r[1] . "\t" . $r[3] . "\t" . $r[2] . "\t" . $r[4] . "\t" . $r[5] . "\t" . $misa{$id} . "\n";
}

exit;
