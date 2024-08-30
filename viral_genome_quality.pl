#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use P3DataAPI;
use JSON::XS;
use File::Slurp;

#
#   Adding in some notes here. The current GTO definition does not allow for an event id
#   To be tagged to something that isn't a feature.  Since overall genome quality is not a feature
#   I cannot tag the event id to it.  In its current state, this program generates an event id and 
#   returns quality but does not link the two.
#
#   Things that we need to add to the GTO
#      1.  viral segment as an aspect of the contig
#          we need:
#          string viral_segment (im currently using replicon_type)
#
#      2.  In feature_quality_measure, I would consider adding:
#          bool     truncated
#          float    expected_min_len (i would prefer in aa, but could be nts)
#          float    expected_max_len (i would prefer in aa, but could be nts)
#          string   exception (free text description of problem)
#          string   feature_quality ("Good" or "Poor") 
#
#      3.  possibly an analysis event Id for the genome quality and contig data
#
#

my($opt, $usage) = describe_options("%c %o",
				    ["ambiguous|a=f"   => "Fraction of ambiguous bases", { default => 0.01 }],
				    ["input|i=s"       => "Input GTO"],
				    ["output|o=s"      => "Output GTO"],
				    ["prefix|p=s"      => "Genome Quality File Prefix", { default => "Viral_Anno" }],
				    ["json|j=s"        => "Full path to the JSON opts file", {default => "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json"}],
				    ["help|h"          => "Show this help message"]);


print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;
my $prefix = $opt->prefix // "Viral_Anno";

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";
$genome_in->{features}->[0] or die "No features in GTO\n"; 

my $json      = decode_json(read_file($opt->json));
$genome_in or die "Error reading json protein feature data";


#
# Create the GTO analysis event
#
my $event = {
    tool_name => "viral_genome_quality",
    execution_time => scalar gettimeofday,
};

my $event_id = $genome_in->add_analysis_event($event);



## We read the GTO to get the family pssm that was used to call the proteins so that
## we can look up which genes are essential.
## also make sure that the GTO has our annotations.
my %pssm_fam;
for my $i (0 .. $#{$genome_in->{features}}) 
{
	if (($genome_in->{features}->[$i]->{type} =~ /CDS/) && ($genome_in->{features}->[$i]->{family_assignments}->[0]->[3] =~ /annotate_by_viral_pssm/))
	{
		$pssm_fam{$genome_in->{features}->[$i]->{family_assignments}->[0]->[0]}++;
	}
}
die "More than one viral family of PSSMs in GTO\n" if scalar(keys %pssm_fam) > 1;
my $fam = (keys %pssm_fam)[0];
$fam or die "GTO has no annotations from annotate_by_viral_pssm tool\n"; 


# Next, we read the annotation json file to initiate a hash of which proteins are essential.
# In this case, essential really means that it is one of the proteins we are looking for.
# All of the proteins are probably essential, but the ones that are variable in length or
# presence/absence in the family are passed over.


my $essential = {};
my %nonessential;  
# I am keeping track of CDs where the segment is declared to annotate segements with 
# and as a double check on segment copy number.

foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	my $anno     = $json->{$fam}->{features}->{$prot}->{anno};

	if (($json->{$fam}->{features}->{$prot}->{copy_num}) && ($json->{$fam}->{features}->{$prot}->{feature_type} =~ /CDS/))
	{
		$essential->{$anno}->{copy_num} = $json->{$fam}->{features}->{$prot}->{copy_num};
		$essential->{$anno}->{segment}  = $json->{$fam}->{features}->{$prot}->{segment};
		$essential->{$anno}->{max_len}  = $json->{$fam}->{features}->{$prot}->{max_len};
		$essential->{$anno}->{min_len}  = $json->{$fam}->{features}->{$prot}->{min_len};
		$essential->{$anno}->{found_copy_num} = 0;
	}
	elsif (($json->{$fam}->{features}->{$prot}->{feature_type} =~ /CDS/i) && ($json->{$fam}->{features}->{$prot}->{segment}) && ($json->{$fam}->{features}->{$prot}->{anno}))
	{
		$nonessential{$anno} = $json->{$fam}->{features}->{$prot}->{segment};
	}
}
unless ($essential){die "No essential data in annotation JSON file\n";}

open (OUT1, ">$prefix.feature_quality");

my $contig_eval = {};
#Next we read the gto and look for the essential proteins
for my $i (0 .. $#{$genome_in->{features}}) 
{
	my $anno = $genome_in->{features}->[$i]->{function};
	my $prot = $genome_in->{features}->[$i]->{protein_translation};
	$prot =~ s/\*$//g; 

	if (exists $essential->{$anno})
	{
		$essential->{$anno}->{found_len} = length ($prot); 
		$essential->{$anno}->{found_contig} -> {$genome_in->{features}->[$i]->{location}->[0]->[0]} ++;
		$essential->{$anno}->{found_copy_num} ++;
		@{$essential->{$anno}->{found_ids}} = $genome_in->{features}->[$i]->{id};
	}
	elsif (exists $nonessential{$anno})  #grab the contig ids for non-essential genes
	{
		$contig_eval->{$nonessential{$anno}}->{$genome_in->{features}->[$i]->{location}->[0]->[0]} ++; 		
	}
}

my $genome_quality = "Good"; 
print OUT1 "ID(s)\tFunction\tProt_Len\tMin_Len\tMax_Len\tExp_Copy_Num\tFound_Copy_Num\tException\n";

foreach (keys %{$essential})
{
	my $exception;
	my $anno = $_; 
	if ($essential->{$anno}->{found_len} > $essential->{$anno}->{max_len}) 
	{
		$exception = "Protein is longer than expected";
	}
	elsif ($essential->{$anno}->{found_len} < $essential->{$anno}->{min_len}) 
	{
		$exception = "Protein is shorter than expected";
	}
	elsif ($essential->{$anno}->{found_copy_num} > $essential->{$anno}->{copy_num}) 
	{
		$exception = "Too many copies of protein";
	}
	elsif ($essential->{$anno}->{found_copy_num} < $essential->{$anno}->{copy_num}) 
	{
		$exception = "Too few copies of protein";
	}
	
	if ($exception){$genome_quality = "Poor";}

	#here, since i am cycling through each essential, 
	#I will populate a hash of the Segment => Contig ID.
	
	foreach (keys %{$essential->{$anno}->{found_contig}})
	{
		my $contig = $_;
		$contig_eval->{$essential->{$anno}->{segment}}->{$contig} ++; 
	}
	
	my $id_string = join (";", @{$essential->{$anno}->{found_ids}});
	print OUT1 "$id_string\t$anno\t$essential->{$anno}->{found_len}\t$essential->{$anno}->{min_len}\t$essential->{$anno}->{max_len}\t$essential->{$anno}->{copy_num}\t$essential->{$anno}->{found_copy_num}\t$exception\n"; 
}
close OUT1;


open (OUT2, ">$prefix.contig_quality");
print OUT2 "Contig(s)\tSegment\tCopy_Num\tLen\tMin_Len\tMax_len\tAMB_Bases\tException\n"; 

#first work through the contigs containing on essential proteins
foreach (keys %{$contig_eval})
{
	my $contig_exception;
	my $seg = $_; 
	my @contigs = keys %{$contig_eval->{$seg}}; 
	my $n_contigs = scalar @contigs;
	my $contig_string = join (";", @contigs); 
	if ($n_contigs > 1)
	{
		$contig_exception = "Too many contigs for $seg";
	}
	if ($n_contigs < 1)# only the essentials are placed into the contig eval hash unpopulated 
	{
		$contig_exception = "Too few contigs for $seg";
	}
	
	my $len;
	my $amb_bases;
	my $frac_amb;
	for my $i (0 .. $#{$genome_in->{contigs}})
	{ 
		if($genome_in->{contigs}->[$i]->{id} eq $contig_string) ### so it will not add the segment if it is broken.
		{
			$len = length($genome_in->{contigs}->[$i]->{dna});
			
			#count the total number of ambiguous bases in the dna string of the contig from the gto.
			$amb_bases = () = $genome_in->{contigs}->[$i]->{dna} =~ /[^atgc]/gi;
			$frac_amb = ($amb_bases/$len);
			
			#add seg and replicon geometry to GTO
			$genome_in->{contigs}->[$i]->{replicon_type} = $seg;

			if ($json->{$fam}->{segments}->{$seg})
			{
				$genome_in->{contigs}->[$i]->{replicon_geometry} = $json->{$fam}->{segments}->{$seg}->{replicon_geometry};
			}
			if (($json->{$fam}->{segments}->{$seg}) && ($len < $json->{$fam}->{segments}->{$seg}->{min_len}))
			{
				$contig_exception = "Segment: $seg, Contig: $contig_string is too short";
			}
			if (($json->{$fam}->{segments}->{$seg}) &&  ($len > $json->{$fam}->{segments}->{$seg}->{max_len}))
			{
				$contig_exception = "Segment: $seg, Contig: $contig_string is too long";
			}		
			if ($frac_amb > ($opt->ambiguous))
			{
				$contig_exception = "Segment: $seg, Contig: $contig_string has more than 1% ambiguous bases";
			}
		}
	}	 
	
	if ($contig_exception){$genome_quality = "Poor";}
	print OUT2 "$contig_string\t$seg\t$n_contigs\t$len\t$json->{$fam}->{segments}->{$seg}->{min_len}\t$json->{$fam}->{segments}->{$seg}->{max_len}\t$amb_bases\t$frac_amb\t$contig_exception\n";
}

print "$genome_in->{id}\t$genome_quality\n";

#shove quality into the GTO, and write the new GTO.
$genome_in->{quality}->{genome_quality} = $genome_quality;
$genome_in->destroy_to_file($opt->output);



