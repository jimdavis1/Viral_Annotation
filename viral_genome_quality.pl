#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use P3DataAPI;
use JSON::XS;
use File::Slurp;


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
				    ["ambiguous|a=f"   => "Fraction of ambiguous bases, (Default = 0.01)", { default => 0.01 }],
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



# Create the GTO analysis event

my $event = {
    tool_name => "viral_genome_quality",
    execution_time => scalar gettimeofday,
};

my $event_id = $genome_in->add_analysis_event($event);



# We read the GTO to get the family pssm that was used to call the proteins so that
# we can look up which genes are essential.
# also make sure that the GTO has our annotations.
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
my %anno_count;

# I am keeping track of CDs where the segment is declared to annotate segements with 
# and as a double check on segment copy number.

foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	my $anno     = $json->{$fam}->{features}->{$prot}->{anno};
	if (($json->{$fam}->{features}->{$prot}->{copy_num}) && ($json->{$fam}->{features}->{$prot}->{feature_type} =~ /(CDS|mat_peptide)/)) 
	{
		$essential->{$anno}->{copy_num} = $json->{$fam}->{features}->{$prot}->{copy_num};
		$essential->{$anno}->{segment}  = $json->{$fam}->{features}->{$prot}->{segment};
		$essential->{$anno}->{max_len}  = $json->{$fam}->{features}->{$prot}->{max_len};
		$essential->{$anno}->{min_len}  = $json->{$fam}->{features}->{$prot}->{min_len};
		$anno_count{$anno} = 0;
	}
	#also grab nonessential proteins with a known segment.
	elsif (($json->{$fam}->{features}->{$prot}->{feature_type} =~ /(CDS|mat_peptide)/i) && ($json->{$fam}->{features}->{$prot}->{segment}) && ($json->{$fam}->{features}->{$prot}->{anno}))
	{
		$nonessential{$anno} = $json->{$fam}->{features}->{$prot}->{segment};
	}
}

unless ($essential){die "No essential data in annotation JSON file\n";}
my $seg_count = {};
my $feature_eval = {};
my $genome_quality = "Good"; 



#Next we read the gto and look for the essential proteins
for my $i (0 .. $#{$genome_in->{features}}) 
{
	if ($genome_in->{features}->[$i]->{type} =~ /(CDS|mat_peptide)/)
	{
		my $anno = $genome_in->{features}->[$i]->{function};
		my $prot = $genome_in->{features}->[$i]->{protein_translation};
		$prot =~ s/\*$//g; 
		my $id = $genome_in->{features}->[$i]->{id};
		my $len = length($prot);

		$feature_eval->{$id}->{ANNO}       = $anno;
		$feature_eval->{$id}->{LEN}        = $len;
		$feature_eval->{$id}->{CONTIG}     = $genome_in->{features}->[$i]->{location}->[0]->[0]; 
		## need to loop through his a second time.
	
	
		#search for exceptions;
		if ($essential->{$anno})
		{
			$feature_eval->{$id}->{MAX_LEN} = $essential->{$anno}->{max_len};
			$feature_eval->{$id}->{MIN_LEN} = $essential->{$anno}->{min_len};
			$feature_eval->{$id}->{EXP_CN}  = $essential->{$anno}->{copy_num};
			$feature_eval->{$id}->{EVALUATED}  = 1;
			$anno_count{$anno} ++;
		
			if ($len > $essential->{$anno}->{max_len}) 
			{
				push @{$feature_eval->{$id}->{EXCEPTION}},   "Protein is longer than expected";
			}
			elsif ($len < $essential->{$anno}->{min_len}) 
			{
				push @{$feature_eval->{$id}->{EXCEPTION}},   "Protein is shorter than expected";
			}
		
			#contigs per segment (essential features)
			$seg_count->{$essential->{$anno}->{segment}}->{$genome_in->{features}->[$i]->{location}->[0]->[0]} ++; 
		}
		#contigs per segment (nonessential features)
		elsif (exists $nonessential{$anno})  #grab the contig ids for non-essential genes
		{
			$seg_count->{$nonessential{$anno}}->{$genome_in->{features}->[$i]->{location}->[0]->[0]} ++; 		
		}
	}
}


# Find Copy number exceptions:

open (OUT1, ">$prefix.feature_quality");
print OUT1 "ID\tFunction\tEvaluated\tProt_Len\tMin_Len\tMax_Len\tCopy_Num\tExp_Copy_Num\tExceptions\n";

foreach (sort keys %{$feature_eval})
{
	my $id     = $_;
	my $anno   = $feature_eval->{$id}->{ANNO};
	my $exp_cn = $feature_eval->{$id}->{EXP_CN};
	my $cn = $anno_count{$anno}; 
	$feature_eval->{$id}->{CN} = $cn;  
	if (($cn > $exp_cn) && ($feature_eval->{$id}->{EVALUATED}))
	{
		push @{$feature_eval->{$id}->{EXCEPTION}}, "Too many copies of protein";
	}
	elsif (($cn < $exp_cn) && ($feature_eval->{$id}->{EVALUATED}))
	{
		push @{$feature_eval->{$id}->{EXCEPTION}}, "Too few copies of protein";
	}
	
	my $exs;
	if (exists $feature_eval->{$id}->{EXCEPTION})
	{
		$exs = join (";", @{$feature_eval->{$id}->{EXCEPTION}});
		$genome_quality = "Poor";
	}
	print OUT1 "$id\t$anno\t$feature_eval->{$id}->{EVALUATED}\t$feature_eval->{$id}->{LEN}\t$feature_eval->{$id}->{MIN_LEN}\t$feature_eval->{$id}->{MAX_LEN}\t$cn\t$exp_cn\t$exs\n";
}

#go back an look for annotations that were missed.
foreach (keys %anno_count)
{	
	if ($anno_count{$_} < 1)
	{
		print OUT1 "\n##\nMissing\t$_\n";
		$genome_quality = "Poor";
	}
}
close OUT1;


#go grap the name of each segement that should be there.
my %segs_needed;
foreach (keys %{$json->{$fam}->{segments}})
{
	$segs_needed{$_} = 0;
}


###Evaluate the Contig Quality

open (OUT2, ">$prefix.contig_quality");
print OUT2 "Contig\tSegment\tCopy_Num\tLen\tMin_Len\tMax_len\tAMB_Bases\tFrac_AMB\tExceptions\n"; 
my $contig_eval = {};



#first work through the contigs containing on essential proteins
foreach (sort keys %{$seg_count})
{
	my $seg = $_; 
	$segs_needed{$seg}++;
	my @contigs = keys %{$seg_count->{$seg}}; 
	my $n_contigs = scalar @contigs;

	foreach (@contigs)
	{
		my $contig = $_;
		$contig_eval->{$contig}->{SEG} = $seg;
			
		#add an exception if count is too high count;
		if ($n_contigs > 1)
		{
			push @{$contig_eval->{$contig}->{EXCEPTION}}, "Too many contigs for $seg";
		}
	
		for my $i (0 .. $#{$genome_in->{contigs}})
		{ 
			if($genome_in->{contigs}->[$i]->{id} eq $contig) 
			{
				
				my $len = length($genome_in->{contigs}->[$i]->{dna});
				$contig_eval->{$contig}->{LEN} = $len; 
			
				#count the total number of ambiguous bases in the dna string of the contig from the gto.
				my $amb_bases = () = $genome_in->{contigs}->[$i]->{dna} =~ /[^atgc]/gi;
				$contig_eval->{$contig}->{NUM_AMB} = ($amb_bases);
				$contig_eval->{$contig}->{FRAC_AMB} = ($amb_bases/$len);
			
				#add seg and replicon geometry to GTO
				$genome_in->{contigs}->[$i]->{replicon_type} = $seg;

				
				if ($json->{$fam}->{segments}->{$seg})
				{
					$genome_in->{contigs}->[$i]->{replicon_geometry} = $json->{$fam}->{segments}->{$seg}->{replicon_geometry};
					$contig_eval->{$contig}->{MIN_LEN} = $json->{$fam}->{segments}->{$seg}->{min_len};
					$contig_eval->{$contig}->{MAX_LEN} = $json->{$fam}->{segments}->{$seg}->{max_len};
				}
				
				## Look for length and quality exceptions:
				if (($json->{$fam}->{segments}->{$seg}) && ($len < $json->{$fam}->{segments}->{$seg}->{min_len}))
				{
					push @{$contig_eval->{$contig}->{EXCEPTION}}, "Contig is too short";				
				}
				if (($json->{$fam}->{segments}->{$seg}) &&  ($len > $json->{$fam}->{segments}->{$seg}->{max_len}))
				{
					push @{$contig_eval->{$contig}->{EXCEPTION}}, "Contig is too long";				
				}		
				if ($contig_eval->{$contig}->{FRAC_AMB} > ($opt->ambiguous))
				{
					push @{$contig_eval->{$contig}->{EXCEPTION}}, "Contig has more than 1% ambiguous bases";
				}

		
			my $exs;
			if (exists $contig_eval->{$contig}->{EXCEPTION})
			{
				$exs = join (";", @{$contig_eval->{$contig}->{EXCEPTION}});
				$genome_quality = "Poor";
			}

			print OUT2 "$contig\t$seg\t$n_contigs\t$len\t$json->{$fam}->{segments}->{$seg}->{min_len}\t$json->{$fam}->{segments}->{$seg}->{max_len}\t$amb_bases\t$contig_eval->{$contig}->{FRAC_AMB}\t$exs\n";

			}
		}	 
	}
}

# create an exception for missing segment.
foreach (keys %segs_needed)
{
	if ($segs_needed{$_} < 1)
	{
		print OUT2 "\n##\nMissing\t$_\n";
		$genome_quality = "Poor";
	}		
}
close OUT2;


print "$genome_in->{id}\t$genome_quality\n";

#shove quality into the GTO, and write the new GTO.
$genome_in->{quality}->{genome_quality} = $genome_quality;
$genome_in->destroy_to_file($opt->output);



