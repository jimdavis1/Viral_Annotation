#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
use JSON::XS;
use gjoseqlib;


my $usage = 'gtos-to-alis-and-trees.pl  [options] -i GTO_DIR

	-h help
	-i GTO DIR
	-e exclude_file list of protein ids with exceptions to exclude.
	-g good quality only [default = considers all gtos]
	-m metadata file name
	-c color file name

	gtos in the directory must all have the file prefix .gto.
	
	
	
';

my ($help, $dir, $exclude_file, $good_only);
my $name_file  = "id.names";
my $color_file = "id.colors";

my $opts = GetOptions( 'h'         => \$help,
                       'i=s'       => \$dir,
                       'g'         => \$good_only,
                       'e=s'       => \$exclude_file,
                       'm=s'       => \$name_file,
                       'c=s'       => \$color_file,
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


open (OUT, ">$name_file"); 
open (CLR, ">$color_file"); 



my $hash = {};

foreach (@files)
{	
	my $file = $_; 
	my $genome_in = GenomeTypeObject->create_from_file("$dir/$file");
	$genome_in or warn "Error reading and parsing input for $file\n";
	my $gid =  $genome_in->{"id"}; 
	my $qual = $genome_in->{"quality"}->{"genome_quality"};
	my $name = $genome_in->{"scientific_name"};
	
	if (($qual =~ /Good/) || (($qual =~ /Poor/) && (! defined $good_only)))
	{
		my @features = @{$genome_in->{"features"}};
		for my $i (0..$#features)
		{
			my $type = $genome_in->{"features"}->[$i]->{"type"};
			my $id   = $genome_in->{"features"}->[$i]->{"id"};
			my $aa   = $genome_in->{"features"}->[$i]->{"protein_translation"};
			my $anno = $genome_in->{"features"}->[$i]->{"function"};
		
			if (($type =~ /CDS/)||($type =~ /mat_peptide/))
			{
				# had to stick this here because family is only recorded for blast-based features.
				my $fam  = $genome_in->{"features"}->[$i]->{"family_assignments"}->[0]->[0];

				unless (exists $exclude{$id})
				{
					$hash->{$anno}->{$id}->{"AA"} = $aa;
					$hash->{$anno}->{$id}->{"NAME"} = $name;
					$hash->{$anno}->{$id}->{"QUAL"} = $qual;
					print OUT "$id\t$name \[$qual\]\n";
					print CLR "$id\t$fam\n";
				}
			}
		}
	}
}

foreach (sort keys %$hash)
{
	my $anno = $_; 
	my $file = $anno;
	$file =~ s/ /_/g;
	$file =~ s/\'//g; 
	$file =~ s/\(//g; 
	$file =~ s/\)//g; 
    open (OUT2, ">$file.aa"); 
	
	my @seqs;
	foreach (keys %{$hash->{$anno}})
	{
		my $id = $_;
		my $name = $hash->{$anno}->{$id}->{"NAME"};
		my $aa = $hash->{$anno}->{$id}->{"AA"};
		my $qual = $hash->{$anno}->{$id}->{"QUAL"};
		$aa =~ s/\*$//g;
		
		my $des = "$anno"."##"."[$name]"."##"."[$qual]";
		push @seqs, ([$id, $des, $aa]);
	}	
	
	&gjoseqlib::print_alignment_as_fasta(\*OUT2, @seqs);
	system "mafft --thread 24 --reorder $file.aa >$file.fa";
	system "FastTree -wag <$file.fa>$file.nwk";
	system "svr_tree_to_html -c $color_file -raw -bar -a $name_file <$file.nwk>$file.html";

}


































