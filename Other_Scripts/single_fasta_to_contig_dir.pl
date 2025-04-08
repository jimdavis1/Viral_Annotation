#! /usr/bin/env perl
use strict;
use Getopt::Long;
use gjoseqlib;


my $usage = 'Single_fasta_to_contig_dir.pl -i Contigs Dir <fasta 

	This program reads a fasta file and makes a contigs directory of separate fasta files for processing.
	
	Options:
		-h help
		-c Contigs Directory name [d = Contigs]

';


my $contigs = "Contigs"; 
my $help;

my $opts = GetOptions( 'h'         => \$help,
                       'c=s'       => \$contigs,
                           
); 


if ($help){die "$usage\n";}
mkdir $contigs;

my @seqs = &gjoseqlib::read_fasta();

for my $i (0..$#seqs)
{
	my $id = $seqs[$i][0];
	open (OUT, ">$contigs/$id.contigs");
	&gjoseqlib::print_alignment_as_fasta(\*OUT, @seqs[$i]); 



}