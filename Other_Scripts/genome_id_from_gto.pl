#! /usr/bin/env perl
use strict;
use GenomeTypeObject;

my $usage = 'genome_id_from_gto.pl genome.gto

    Reads a GTO file and prints the filename and genome ID tab delimited.

';

my $file = $ARGV[0] or die $usage;

my $genome = GenomeTypeObject->create_from_file($file);
$genome or die "Error reading and parsing GTO from $file\n";

my $gid = $genome->{"id"};

print "$file\t$gid\n";