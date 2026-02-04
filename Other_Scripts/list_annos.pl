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
	print "| Taxon | PSSM_Name | Anno | Feat_Type | Copy_Number | Segment | Genome_Quality |\n";
	print "|:------|:---------:|:----:|:---------:|:-----------:|:-------:|:--------------:|\n";	
}
else
{
	print "Taxon\tPSSM_Name\tAnno\tFeat_Type\tCopy_Number\tSegment\tGenome_Quality\n";
}


foreach (sort keys %$json)
{
	my $fam = $_;
	
	foreach (sort keys %{$json->{$fam}->{features}})
	{
		my $pssm = $_;
		my $feat = $json->{$fam}->{features}->{$pssm};

		my $type = $feat->{feature_type};

		# anno: prefer anno, fall back to new_anno
		my $anno;
		if ($feat->{anno})
		{
			$anno = $feat->{anno};
		}
		elsif ($feat->{new_anno})
		{
			$anno = $feat->{new_anno};
		}

		# copy_num: use value if present, otherwise empty string
		my $copy_num = defined $feat->{copy_num} ? $feat->{copy_num} : '';

		# segment: use value if present, otherwise empty string
		my $segment = defined $feat->{segment} ? $feat->{segment} : '';

		# genome_quality: 1 if BOTH min_len and max_len are defined, 0 otherwise
		my $genome_quality = (defined $feat->{min_len} && defined $feat->{max_len}) ? 1 : 0;

		if ($markdown)
		{
			print "| $fam | $pssm | $anno | $type | $copy_num | $segment | $genome_quality |\n";
		}
		else
		{
			print "$fam\t$pssm\t$anno\t$type\t$copy_num\t$segment\t$genome_quality\n";	
		}
	
	}
}