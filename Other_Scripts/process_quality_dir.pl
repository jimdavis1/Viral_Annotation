#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
#use P3DataAPI;
use JSON::XS;


my $usage = 'process_quality_dir.pl -i Quality_Dir -p file_prefix

		-h   help
		-i   GTO quality DIR
		-p   File prefix for outputs [default = "Output"]
			 Prefix.qual ([Tuple of original id, new id, and quality])
			 Prefix.stats
	    -c   Contigs file directory.  This is used to determine the contigs
	         that were downloaded but not processed due to blast dissimilarity

	The quality directory must be formatted with three files per genome

		1.  genome_id.qual.gto
		2.  genome_id.contig_quality
		3.  genome_id.feature_quality
	
	Genome counts and fractions are based on the total number of gtos in the directory.

';

my ($help, $dir, $contigsD);
my $prefix = "Output"; 
my $opts = GetOptions( 'h'         => \$help,
                       'i=s'       => \$dir,
                       'c=s'       => \$contigsD,
                       'p=s'       => \$prefix); 


if ($help){die "$usage\n";}
unless ($dir){die "Must declare a directory with genome quality data\n"; }

opendir (DIR, "$dir");
my @files = grep{$_ !~ /^\./}readdir(DIR); 
close DIR;

### Get the diverse genomes that weren't processed due to blast dissimilarity

my @cons;
my @diverse;
if ($contigsD)
{
	opendir (DIR, "$contigsD");
	open (OUT, ">$prefix.diverse.no_analysis");

	my @cons = grep{$_ !~ /^\./}readdir(DIR); 
	close DIR;
	my %processed; 
	foreach (@files)
	{	
		my $qual_id = $_; 
		$qual_id =~ s/\.qual.gto//g; 
		$processed{$qual_id} = 1;
	}
	foreach (@cons)
	{
		my $genome = $_;
		$genome =~ s/\.contig.+//g; 
		unless (exists $processed{$genome})
		{
			push @diverse, $genome; 
			print OUT "$genome\n"; 
		}
	}	
}



#
# First, I am going to read the GTO and get good vs. poor
#
my @good; 
my @poor; 
my @bad_contigs;
my @bad_features;
my %bad_genomes_by_feature;
my %contig_exceptions;
my $feature_exceptions = {};
my %bad_contig_files;
foreach (@files)
{
	if ($_ =~ /\.gto$/)
	{
		print STDERR "$_\n"; 
		my $file = "$dir/$_";
		my $genome_in = GenomeTypeObject->create_from_file($file);
		$genome_in or warn "$file:  Error reading and parsing input\n";
		my $qual = $genome_in->{quality}->{genome_quality};
		my $orig_id = $_;
		$orig_id =~ s/\.gto//g; 
		my $new_id = $genome_in->{id};
		if ($qual =~ /Good/)
		{
			push @good, ([$orig_id, $new_id, $qual]);
		}
		elsif ($qual =~ /Poor/)
		{
			push @poor, ([$orig_id, $new_id, $qual]);
		}		
	}
	elsif ($_ =~ /\.contig_quality$/)
	{
		open (IN, "<$dir/$_"), or warn "could not open $dir/$_\n"; 
		{
			my $file = $_; 
			$file =~ s/\.contig_quality//g; 
			while (<IN>)
			{
				chomp;
				unless ($_ =~ /^Contig/)
				{
					my ($contig, $segment, $copynum, $len, $min, $max, $amb, $fracamb, $exc) = split /\t/; 
					if ($exc)
					{
						push @bad_contigs, ([$file, $contig, $segment, $copynum, $len, $min, $max, $amb, $fracamb, $exc]);
						$contig_exceptions{$exc}++;
						$bad_contig_files{$file} ++;
					}
				}
			}
		}
		close IN;
	}
	elsif ($_ =~ /\.feature_quality$/)
	{
		open (IN, "<$dir/$_"), or warn "could not open $dir/$_\n"; 
		{
			my $file = $_; 
			$file =~ s/\.feature_quality//g; 
			while (<IN>)
			{
				chomp;
				unless ($_ =~ /^ID/)
				{
					my ($id, $fun, $eval, $len, $min, $max, $cn, $ecn, $exc) = split /\t/; 
					if ($exc)
					{
						$bad_genomes_by_feature{$file} ++;
						
						unless (exists $bad_contig_files{$file})
						{
							push @bad_features, ([$file, $id, $fun, $eval, $len, $min, $max, $cn, $ecn, $exc]);
							$feature_exceptions->{$fun}->{$exc}++;
						}
					}
				}
			}
		}
	}
}












open (OUT, ">$prefix.qual");
for my $i (0..$#poor)
{	
	print OUT "$poor[$i][0]\t$poor[$i][1]\t$poor[$i][2]\n"; 
}
for my $i (0..$#good)
{	
	print OUT "$good[$i][0]\t$good[$i][1]\t$good[$i][2]\n"; 
}
close OUT; 

open (OUT, ">$prefix.bad_contigs");
for my $i (0..$#bad_contigs)
{	
	my $arrayR= $bad_contigs[$i]; 
	my @array = @$arrayR;
	print OUT join "\t", @array, "\n";  
}
close OUT;


my $sum_good = scalar @good;
my $sum_poor = scalar @poor;
my $sum = ($sum_good + $sum_poor); 

my $frac_good = 0;
my $frac_poor = 0;

my $sum_bad_contigs = scalar @bad_contigs;
my $frac_bad_contigs = 0;

my $sum_bad_by_feature = keys %bad_genomes_by_feature;
my $frac_bad_by_feature = 0;

if ($sum_good)
{
	$frac_good = sprintf("%.3f", ($sum_good/$sum));
}
if ($sum_poor)
{
	$frac_poor = sprintf("%.3f", ($sum_poor/$sum));
}
if ($sum_bad_contigs)
{
	$frac_bad_contigs = sprintf("%.3f", ($sum_bad_contigs/$sum));
}
if ($sum_bad_by_feature)
{
	$frac_bad_by_feature = sprintf("%.3f", ($sum_bad_by_feature/$sum));
}



print STDERR "Total = $sum\tGood = $sum_good \($frac_good\)\tPoor = $sum_poor \($frac_poor\)\n";
print STDERR "Bad Genomes by contigs = $sum_bad_contigs \($frac_bad_contigs\)\n"; 
print STDERR "Bad Genomes by feature (all genomes) = $sum_bad_by_feature \($frac_bad_by_feature\)\n"; 










print STDERR "\nContig Exceptions:\n"; 
foreach (sort{$contig_exceptions{$b} <=> $contig_exceptions{$a}}keys %contig_exceptions)
{
	print STDERR "\t$contig_exceptions{$_}\t$_\n"; 
}

print STDERR "\nFeature Exceptions (not counting genomes with contig exceptions):\n"; 
foreach (keys %$feature_exceptions)
{
	my $feat = $_;
	foreach (keys %{$feature_exceptions->{$feat}})
	{
		print STDERR "$feat\t$_\t$feature_exceptions->{$feat}->{$_}\n";
	}
}



























