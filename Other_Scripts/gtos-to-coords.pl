#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
use JSON::XS;
use gjoseqlib;


my $usage = 'gtos-to-coords.pl  [options] -i GTO_DIR

	-h help
	-i GTO DIR
	-e exclude_file list of filenames to exclude

	gtos in the directory must all have the file prefix .gto.
	This purges all poor-quality genomes.
	
	
';

my ($help, $dir, $exclude_file);


my $opts = GetOptions( 'h'         => \$help,
                       'i=s'       => \$dir,
                       'e=s'       => \$exclude_file,
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


foreach (@files)
{	
	my $file = $_;
	my $old_id = $file;
	$old_id =~ s/\.qual.+//g; 
	
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
			my $type   = $genome_in->{"features"}->[$i]->{"type"};		
			my $peg_id = $genome_in->{"features"}->[$i]->{"id"};		
			my $anno   = $genome_in->{"features"}->[$i]->{"function"};	
			my @locs   = @{$genome_in->{"features"}->[$i]->{"location"}};	
			if (scalar @locs > 1)  # no broken features.
			{
				warn "$file:  $peg_id:  multi-location feature";
				next;
			}
			else
			{
				my $contig = $genome_in->{"features"}->[$i]->{"location"}->[0]->[0];
				my $begin  = $genome_in->{"features"}->[$i]->{"location"}->[0]->[1];
				my $strand = $genome_in->{"features"}->[$i]->{"location"}->[0]->[2];
				my $len    = $genome_in->{"features"}->[$i]->{"location"}->[0]->[3];
				my $end; 
				if ($strand =~ /\+/)
				{
					$end = (($begin + $len ) - 1);
				}
				elsif ($strand =~ /\-/)
				{
					$end = (($begin - $len) +1);
				}
			
				print "$old_id\t$gid\t$peg_id\t$begin\t$end\t$anno\t$name\n";
			}
		}
	}
}



































