#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
use JSON::XS;
use gjoseqlib;


my $usage = 'gtos-dir-to-wgs-tree-single-seg.pl  [options] -i GTO_DIR

	-h help
	-i GTO DIR
	-e exclude_file list of filenames to exclude
	-p prefix [D = Viral_WGS]

	gtos in the directory must all have the file prefix .gto.
	This purges all poor-quality genomes.
	
	
';

my ($help, $dir, $exclude_file);
my $prefix = "Viral_WGS";


my $opts = GetOptions( 'h'         => \$help,
                       'i=s'       => \$dir,
                       'e=s'       => \$exclude_file,
                       'p=s'       => \$prefix
                       ); 


if ($help){die "$usage\n";}
unless ($dir){die "Must declare a directory with GTOs"; }

# set up the exclude hash.
my %exclude;
if ($exclude_file)
{
	open (IN, "<$exclude_file"), or die "cannot open exclude file\n";
	%exclude = map{chomp; $_, 0}(<IN>);
	close IN;
}

# read in the GTO dir.
opendir (DIR, "$dir");
my @files = grep{$_ =~ /\.gto/}readdir(DIR); 
close DIR;


open (CLR, ">$prefix.colors"); 
open (NAME, ">$prefix.name"); 

my $hash = {};
my @seqs;

foreach (@files)
{	
	my $file = $_; 
	unless (exists $exclude{$file})
	{
	
		my $genome_in = GenomeTypeObject->create_from_file("$dir/$file");
		$genome_in or warn "Error reading and parsing input for $file\n";
		my $gid =  $genome_in->{"id"}; 
		my $name = $genome_in->{"scientific_name"};
		my $qual = $genome_in->{"quality"}->{"genome_quality"}; # i chose to ditch this.
		next if ($qual =~ /Poor/);	
		
		#get the annotated family for the color file. 
		my $fam;
		my @features = @{$genome_in->{"features"}};
		for my $i (0..$#features)
		{
			my $type = $genome_in->{"features"}->[$i]->{"type"};		
			if ($type =~ /CDS/)
			{
				$fam  = $genome_in->{"features"}->[$i]->{"family_assignments"}->[0]->[0];
			}
		}
		my @contigs = @{$genome_in->{"contigs"}};
		if (scalar @contigs > 1)
		{
			warn "$file is multi-segmented\n";
		}
		else
		{
			print CLR "$gid\t$fam\n";
			my $dna = $genome_in->{"contigs"}->[0]->{"dna"};
			my $des = "$name [$fam]";
			print NAME "$gid\t$des\n";
			push @seqs, ([$gid, $des, $dna]);
		}
	}
}
close CLR;
close NAME;

open (OUT, ">$prefix.dna");
&gjoseqlib::print_alignment_as_fasta(\*OUT, @seqs);
system "mafft --thread 24 --reorder $prefix.dna >$prefix.fa";
system "FastTree -nt -gtr <$prefix.fa>$prefix.nwk";
system "svr_tree_to_html -c $prefix.colors -raw -bar -a $$prefix.name <$prefix.nwk>$prefix.html";




































