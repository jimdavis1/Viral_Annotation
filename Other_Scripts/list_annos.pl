#! /usr/bin/env perl
use strict;
use Getopt::Long;
use JSON::XS;
use File::Slurp;
use Data::Dumper;

my $usage = 'list_annos_from_pssms.pl -b base_directory
		
		        -j   Full path to the options file in JSON format which carries data for a match 
		             (D = /home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json)
				-m print markdown instead of tab-delimited
				-h help';

my ($help, $markdown);
my $json_file = "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json"; 		
my $opts = GetOptions( 'h'         => \$help,
                       'j=s'       => \$json_file,
                       'm'         => \$markdown);

if ($help){die "$usage\n";}

open (IN, "<$json_file"), or die "Cant find JSON options file\n"; 
my $json = decode_json(scalar read_file(\*IN));
close IN;

if ($markdown)
{
	print "| Taxon | PSSM_Name | Feat_Type | Anno |\n";
	print "|:------|:---------:|:---------:|:---:|\n";	
}
else
{
	print "Taxon\tPSSM_Name\tFeat_Type\tAnno\n";
}


foreach (sort keys %$json)
{
	my $fam = $_;
	
	foreach (sort keys %{$json->{$fam}->{features}})
	{
		my $pssm = $_;
		my $type = $json->{$fam}->{features}->{$pssm}->{feature_type};
		
		my $anno;
		if ($json->{$fam}->{features}->{$pssm}->{anno})
		{
			$anno = $json->{$fam}->{features}->{$pssm}->{anno};
		}
		elsif ($json->{$fam}->{features}->{$pssm}->{new_anno})
		{
			$anno = $json->{$fam}->{features}->{$pssm}->{new_anno};
		}
	
		if ($markdown)
		{
			print "| $fam | $pssm | $type | $anno |\n";
		}
		else
		{
			print "$fam\t$pssm\t$type\t$anno\n";	
		}
	
	}
}
