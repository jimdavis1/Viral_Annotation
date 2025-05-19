#! /usr/bin/env perl
use strict;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;
use Cwd;
use gjoseqlib;
use Proc::ParallelLoop;
use GenomeTypeObject;


my $usage = 'annotate_viral_taxon.pl [options] -i Taxon_Name

		-h   help
		-i   taxon name
		-t   threads (D = N)
		
		-a   Annotation GTO folder (D = Anno_GTOs)
		-q   Quality GTO folder (D = Quality_GTOs)
		-c   Contig GTO folder (D = Contig_GTO)
		-e   Transcript Edit GTO folder (D = TE_GTO)
        
        -j   Full path to the options file in JSON format 
             (D = /home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json)

		
		-exclude:  File of genome IDs to exclude from the analysis
		-min:  force a min length rather than looking it up in the json
		-max:  force a max length rather than looking it up in the json



';


my $contig_gto   = "Contig_GTOs";
my $anno_gto     = "Anno_GTOs"; 
my $te_gto       = "TE_GTOs";
my $quality_gto  = "Quality_GTOs";
my $contigD      = "Contigs";
my $threads      = 24;
my $json         = "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json"; 
my $base = getcwd;

my ($help, $taxon, $min, $max, $exF);

my $opts = GetOptions( 'h'         => \$help,
                       'i=s'       => \$taxon,
                       't=i'       => \$threads,
                       'j=s'       => \$json,
                       'a=s'       => \$anno_gto,
                       'q=s'       => \$quality_gto,
                       'c=s'       => \$contig_gto,
                       'd=s'       => \$contigD,
                       'e=s'       => \$te_gto,
                       'min=i'     => \$min,
                       'max=i'     => \$max,
                       'exclude=s' => \$exF
                       
                       
                       
); 

if ($help){die "$usage\n";}
unless ($taxon){die "Must declare a taxon name\n"; }
if (($min && !$max)||($max && !$min)){die "must declare both min and max len values as intergers\n";}

my %exclude;
if ($exF)
{
	open (IN, "<$exF"), or die "could not open exclude file\n"; 
	while (<IN>)
	{
		chomp;
		$exclude{$_} ++;
	}
	close IN;
}


open (IN, "<$json"), or die "Cant find JSON options file, use -opt\n"; 
my $options = decode_json(scalar read_file(\*IN));
close IN;

unless ($min && $max)
{
	## Step 1, load the json and get the max and min genome lengths for the taxon
	unless (exists $options->{$taxon}->{segments}){die "No segments for $taxon\n"; }
	
	my $max = 0;
	my $min = 0;
	
	foreach my $segment (values %{$options->{$taxon}->{"segments"}}) 
	{
		$max += $segment->{"max_len"};
		$min += $segment->{"min_len"};
	}
}






## Step 2 Create the Directories
#my $transcript_edit = 0;
#foreach (keys %{$options->{$taxon}->{"features"}})
#{
#	if ($options->{$taxon}->{"features"}->{$_}->{"special"} eq "transcript_edit")
#	{
#		$transcript_edit = 1;
#	}
#}
mkdir $contig_gto;
mkdir $anno_gto;
mkdir $quality_gto;
mkdir $contigD;
mkdir $te_gto;                   


## Step 3
#Gather genomes that are an appropriate length from BVBRC
#This won't work for segmented viruses unless I add a bit of logic linking segments here.
#(which I haven't done)

open (IN, "echo $taxon | query_PATRIC_bob.pl -c genome -i taxon_lineage_names -r \"genome_id genome_length genome_name\" | ");

open (GOOD, ">$taxon.processed"); 
open (BAD, ">$taxon.rejected"); 


my %good_genomes;
while (<IN>)
{
	chomp;
	my ($id, $len, $name) = split /\t/;  	
	
	print "$_\n"; 
	
	unless (exists $exclude{$id})
	{
		if (($len < $max) && ($len > $min))
		{
			$good_genomes{$id} = $name;
			print GOOD "$id\t$len\t$name\n";
		}
		else
		{
			print BAD "$id\t$len\t$name\n";
		}
	}
}

my @good = keys %good_genomes;

close IN;






## Step 4 Download the contig files, create GTO, Annotate, TE, Quality




pareach(
	\@good,  
	
	sub 
	{ 
		my $id = shift @_;
		my $tax = $id;
		$tax =~ s/\..+//g;
		my $name = $good_genomes{$id}; 
		 
		system "echo $id | BVBRC_clean_contigs.pl -d $contigD"; 
		system "rast-create-genome --ncbi-taxonomy-id $tax --scientific-name \"$name\" --domain Viruses --genetic-code 11 --contigs $contigD/$id.contigs >$contig_gto/$id.contig.gto";
		
		chdir $anno_gto;
		
		system "annotate_by_viral_pssm-GTO.pl --threads 4 --prefix $id <$base/$contig_gto/$id.contig.gto > $id.anno.gto";
		chdir $base;
		
		
		#Get the family off of the annotated GTO
		my $genome_in = GenomeTypeObject->create_from_file("$anno_gto/$id.anno.gto");
		$genome_in or warn "Error reading and parsing input for $id\n";
		$genome_in->{features}->[0] or warn "No features in GTO for $id\n"; 

		my %pssm_fam;
		for my $i (0 .. $#{$genome_in->{features}}) 
		{
			if (($genome_in->{features}->[$i]->{type} =~ /CDS/) && ($genome_in->{features}->[$i]->{family_assignments}->[0]->[3] =~ /LowVan Annotate/))
			{
				$pssm_fam{$genome_in->{features}->[$i]->{family_assignments}->[0]->[0]}++;
			}
		}
		warn "More than one viral family of PSSMs in GTO\n" if scalar(keys %pssm_fam) > 1;
		my $fam = (keys %pssm_fam)[0];
		$fam or warn "GTO has no annotations from LowVan Annotate tool\n"; 
		
		my $transcript_edit = 0;
		#read the json opts file to see if there is transcript editing 
		foreach (keys %{$options->{$fam}->{"features"}})
		{
			if ($options->{$fam}->{"features"}->{$_}->{"special"} eq "transcript_edit")
			{
				$transcript_edit = 1;
			}
		}
		
		if ($transcript_edit)
		{
			chdir $te_gto;
			system "get_transcript_edited_features.pl --threads 4 -i $base/$anno_gto/$id.anno.gto -o $id.te.gto";
			chdir $base;
			chdir $quality_gto;
			system "viral_genome_quality.pl -p $id -i $base/$te_gto/$id.te.gto -o $id.qual.gto";		
			chdir $base;
		}
	
		else
		{
  			chdir $quality_gto;
			system "viral_genome_quality.pl -p $id -i $base/$anno_gto/$id.anno.gto -o $id.qual.gto";
			chdir $base;
		}
	}, 
	{ Max_Workers => $threads });






























