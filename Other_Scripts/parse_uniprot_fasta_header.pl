#! /usr/bin/env perl
use strict;
use Getopt::Long;
use gjoseqlib;

my $usage = 'parse_uniprot_fasta_header.pl <UniProt.fasta >Output.fasta
        -h   help
        -t   table output (default = fasta)
        	 [Table format returns FullID\tSwissPRotID\tUniProtID\tANNO\tORG\tTAX\tGeneName\SeqVariant]
             
        
        -sp  give the swissprot or trembl id instead of uniprot (default) in Fasta header 
';

my ($help, $table, $sp_id);
GetOptions(
    'h'    => \$help,
    'sp'   => \$sp_id,
    't'    => \$table
) or die "$usage\n";

if ($help){die "$usage\n";}

my @seqs = &gjoseqlib::read_fasta();

for my $i (0..$#seqs)
{
	my ($swiss, $uni, $full);
	if ($seqs[$i][0] =~ /^(\w+)\|(\w+)\|(\w+)/) 
	{
		$full = $seqs[$i][0];
		($swiss, $uni) = ($2, $3);
	} 
	else 
	{
    	warn "Header string does not match expected format: $seqs[$i][0]\n";
	}

	my $ann = $seqs[$i][1];
	my ($ox, $os, $gn, $pe, $sv);

	$ann =~ s/\w\w\=.+//g;
	if ($seqs[$i][1] =~ /OS\=/)
	{
		$os  = $seqs[$i][1];
		$os  =~ s/.+OS\=//g; 
		$os  =~ s/\s\w\w\=.+//g;
	}
	if ($seqs[$i][1] =~ /OX\=/)
	{
		$ox  = $seqs[$i][1];
		$ox  =~ s/.+OX\=//g; 
		$ox  =~ s/\s\w\w\=.+//g;
	}
	if ($seqs[$i][1] =~ /GN\=/)
	{
		$gn  = $seqs[$i][1];
		$gn  =~ s/.+GN\=//g; 
		$gn  =~ s/\s\w\w\=.+//g;
	}
	if ($seqs[$i][1] =~ /SV\=/)
	{
		$sv  = $seqs[$i][1];
		$sv  =~ s/.+SV\=//g; 
		$sv  =~ s/\s\w\w\=.+//g;
	}

	if ($table)
	{
		print "$full\t$swiss\t$uni\t$ann\t$os\t$ox\t$gn\t$sv\n"; 

	}
	elsif (! $table)
	{
		my $id = $uni;
		if ($sp_id){$id = $swiss};
		my $des = "$ann   [$os]";
		&gjoseqlib::print_alignment_as_fasta([$id, $des, $seqs[$i][2]]);
	}
}

























