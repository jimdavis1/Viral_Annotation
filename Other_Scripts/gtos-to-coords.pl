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
		# skip empty files
		unless (-s "$dir/$file")
		{
			warn "Skipping empty file: $file\n";
			next;
		}

		my $genome_in = eval { GenomeTypeObject->create_from_file("$dir/$file") };
		if ($@ || !$genome_in)
		{
			warn "Error reading and parsing input for $file: " . ($@ || "unknown error") . "\n";
			next;
		}

		my $gid  = eval { $genome_in->{"id"} };
		my $name = eval { $genome_in->{"scientific_name"} };
		my $qual = eval { $genome_in->{"quality"}->{"genome_quality"} };

		unless (defined $gid && defined $name)
		{
			warn "Skipping $file: missing required fields (id or scientific_name)\n";
			next;
		}

		next if (defined $qual && $qual =~ /Poor/);

		#get the annotated family for the color file. 
		my $fam;
		my @features = eval { @{$genome_in->{"features"}} };
		if ($@)
		{
			warn "Skipping $file: could not read features: $@\n";
			next;
		}

		for my $i (0..$#features)
		{
			my ($type, $peg_id, $anno, @locs);
			eval
			{
				$type   = $genome_in->{"features"}->[$i]->{"type"};
				$peg_id = $genome_in->{"features"}->[$i]->{"id"};
				$anno   = $genome_in->{"features"}->[$i]->{"function"};
				@locs   = @{$genome_in->{"features"}->[$i]->{"location"}};
			};
			if ($@)
			{
				warn "Skipping feature $i in $file: $@\n";
				next;
			}

			if (scalar @locs > 1)  # no broken features.
			{
				warn "$file:  $peg_id:  multi-location feature";
				next;
			}
			else
			{
				my ($contig, $begin, $strand, $len);
				eval
				{
					$contig = $genome_in->{"features"}->[$i]->{"location"}->[0]->[0];
					$begin  = $genome_in->{"features"}->[$i]->{"location"}->[0]->[1];
					$strand = $genome_in->{"features"}->[$i]->{"location"}->[0]->[2];
					$len    = $genome_in->{"features"}->[$i]->{"location"}->[0]->[3];
				};
				if ($@ || !defined $begin || !defined $strand || !defined $len)
				{
					warn "Skipping feature $peg_id in $file: bad location data\n";
					next;
				}

				my $end;
				if ($strand =~ /\+/)
				{
					$end = (($begin + $len ) - 1);
				}
				elsif ($strand =~ /\-/)
				{
					$end = (($begin - $len) +1);
				}
				else
				{
					warn "Skipping feature $peg_id in $file: unrecognised strand '$strand'\n";
					next;
				}

				print "$old_id\t$gid\t$peg_id\t$begin\t$end\t$anno\t$name\n";
			}
		}
	}
}