#! /usr/bin/env perl
use strict;
use gjoseqlib;
use Getopt::Long;


my $usage = 'BVBRC_clean_contigs.pl [options] <List of genome ids >ID.contigs

		-d   output directory (defaults to cwd)
		-h   help

';

my ($dir, $help); 

my $opts = GetOptions( 'h'         => \$help,
                       'd=s'       => \$dir);

if ($help){die "$usage\n";}
if ((defined $dir) && (! -d $dir))
{
	mkdir $dir;
}



while (<>)
{
	chomp; 
	my @seqs; 
	my $id = $_; 
	print STDERR "Running: $id\n";

	open (IN, "echo $id | query_PATRIC_bob.pl -c genome_sequence -i genome_id -r \"genome_id accession genome_name sequence\"  |");  
	while (<IN>)
	{
		chomp; 
		my ($id, $acc, $name, $contig) = split /\t/;
		my $des = "$acc $name"; 
		if ((defined $id) && (defined $contig)) 
		{ 

			if ($dir)
			{
				open (OUT, ">$dir/$id.contigs");  
			}
			else
			{
				open (OUT, ">$id.contigs");
			}

			push @seqs, ([$id, $des, $contig]);
		 	&gjoseqlib::print_alignment_as_fasta(\*OUT, @seqs);
		 }
	}
}

	
	 
