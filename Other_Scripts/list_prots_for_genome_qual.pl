#! /usr/bin/env perl
use strict;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;


my $usage = 'list_prots_for_genome_qual.pl [options] Viral-PSSM.json

	This program just reads and parses the Viral-PSSM.json file so that
	We can see which proteins and segments are used for each taxon.
		-h   help
';

my $help;


my $opts = GetOptions( 'h'         => \$help); 

if ($help){die "$usage\n";}

my $file = shift @ARGV;
open (IN, "<$file"), or die "Cant find JSON options file\n"; 
my $json = decode_json(scalar read_file(\*IN));
close IN;

print "Taxon\tAnno\tfeature type\tcopy number\tsegment\tis required\n";
foreach (sort keys %$json)
{
	my $taxon = $_; 
	foreach (sort keys %{$json->{$taxon}->{'features'}})
	{
		my $feat = $_;
		my $anno = $json->{$taxon}->{'features'}->{$feat}->{'anno'};
		my $cn = $json->{$taxon}->{'features'}->{$feat}->{'copy_num'};
		my $seg = $json->{$taxon}->{'features'}->{$feat}->{'segment'};
		my $ft = $json->{$taxon}->{'features'}->{$feat}->{'feature_type'};
		if ($seg && $cn)
		{
			print "$taxon\t$anno\t$ft\t$cn\t$seg\tyes\n";
		}
		elsif (($seg) && (! $cn))
		{
			print "$taxon\t$anno\t$ft\t$cn\t$seg\tno\n";
		}
	}
}