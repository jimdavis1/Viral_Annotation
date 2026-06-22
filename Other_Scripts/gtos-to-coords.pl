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
	-e exclude_file  list of filenames to exclude
	-l list_file     list of genome IDs to process (first column used;
	                 matches ID.qual.gto or ID.gto within GTO_DIR)
	-a               all genomes (disables poor-quality filtering)

	GTOs in the directory must have the file suffix .gto.
	This purges all poor-quality genomes unless -a is set.

';

my ($help, $dir, $exclude_file, $list_file, $all);


my $opts = GetOptions( 'h'   => \$help,
                       'i=s' => \$dir,
                       'e=s' => \$exclude_file,
                       'l=s' => \$list_file,
                       'a'   => \$all,
                     );


if ($help) { die "$usage\n"; }
unless ($dir) { die "Must declare a directory with GTOs\n"; }

# set up the exclude hash
my %exclude;
if ($exclude_file)
{
	open(IN, "<$exclude_file") or die "Cannot open exclude file $exclude_file\n";
	%exclude = map { chomp; $_, 0 } (<IN>);
	close IN;
}

# build the list of files to process
my @files;

if ($list_file)
{
	# read genome IDs from list file (first column only)
	open(IN, "<$list_file") or die "Cannot open list file $list_file\n";
	while (<IN>)
	{
		chomp;
		next if /^\s*$/;
		my $id = (split /\t/, $_)[0];
		$id =~ s/\s+$//;

		# check for ambiguous files (both ID.qual.gto and ID.gto present)
		my @candidates = grep { -e "$dir/$_" } ("$id.qual.gto", "$id.gto");

		if (scalar @candidates == 2)
		{
			warn "Skipping genome ID $id: both $id.qual.gto and $id.gto exist in $dir\n";
			next;
		}
		elsif (scalar @candidates == 1)
		{
			push @files, $candidates[0];
		}
		else
		{
			warn "Could not find GTO file for genome ID $id in $dir (tried $id.qual.gto and $id.gto)\n";
		}
	}
	close IN;
}
else
{
	# default: process all .gto files in the directory
	opendir(DIR, "$dir") or die "Cannot open directory $dir\n";
	my @all_gto_files = grep { /\.gto$/ } readdir(DIR);
	closedir DIR;

	# check for ambiguous IDs (both ID.gto and ID.qual.gto present)
	my %id_to_files;
	for my $f (@all_gto_files)
	{
		my $stripped = $f;
		$stripped =~ s/\.qual\.gto$//;
		$stripped =~ s/\.gto$//;
		push @{ $id_to_files{$stripped} }, $f;
	}

	for my $stripped_id (sort keys %id_to_files)
	{
		if (scalar @{ $id_to_files{$stripped_id} } > 1)
		{
			warn "Skipping genome ID $stripped_id: ambiguous files found: "
			     . join(", ", @{ $id_to_files{$stripped_id} }) . "\n";
		}
		else
		{
			push @files, $id_to_files{$stripped_id}->[0];
		}
	}
}


foreach (@files)
{
	my $file   = $_;
	my $old_id = $file;
	$old_id =~ s/\.qual\.gto$//;
	$old_id =~ s/\.gto$//;

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

		next if (!$all && defined $qual && $qual =~ /Poor/);

		my @features = eval { @{ $genome_in->{"features"} } };
		if ($@)
		{
			warn "Skipping $file: could not read features: $@\n";
			next;
		}

		for my $i (0 .. $#features)
		{
			my ($type, $peg_id, $anno, @locs);
			eval
			{
				$type   = $genome_in->{"features"}->[$i]->{"type"};
				$peg_id = $genome_in->{"features"}->[$i]->{"id"};
				$anno   = $genome_in->{"features"}->[$i]->{"function"};
				@locs   = @{ $genome_in->{"features"}->[$i]->{"location"} };
			};
			if ($@)
			{
				warn "Skipping feature $i in $file: $@\n";
				next;
			}

			if (scalar @locs > 1)    # no broken features
			{
				warn "$file:  $peg_id:  multi-location feature\n";
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
					$end = ($begin + $len) - 1;
				}
				elsif ($strand =~ /\-/)
				{
					$end = ($begin - $len) + 1;
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
