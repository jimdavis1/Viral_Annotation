#! /usr/bin/env perl
use strict; 
use Data::Dumper;
use Getopt::Long;
use gjoseqlib;

my $usage = 'prots_to_alis_and_trees.pl -m metadata results_dir 
		
		This program reads a directory of .faa files from the annotation script
		And turns them into alignments and trees. 
		
						
		-h help
		-m metadata_file
		-p prefix (name to append to each output file)
		-x remove sequences with ambiguous residues 

		
		This program requires tab-delimited metadata a file in the following format:
		1. new genome id 
		2. genome name
		3. family
		4. genome id made by the annotation script
		5. blank
		6. contig id
		7. source
		6. feature type
		8. feature id
		9. start
	   10. stop
	   11. strand
	   12. pssm
	   13. anno
 		
 	   # I generate this file by appending the first 3 columns to the feature table
 	     that is currently produced by the annotation script.	 		 		

';

my ($help, $metadata, $prefix, $remove_x);
my $opts = GetOptions( 'h'         => \$help,
                       'm=s'       => \$metadata,
                       'p=s'       => \$prefix,
                       'x'         => \$remove_x); 

if ($help){die "$usage\n";}
if ((! $metadata) || (! $prefix)){die "must declare -p prefix and -m metadata file\n";}
my $anno_dir = shift @ARGV;

opendir (DIR, "$anno_dir");
my @files = grep{$_ =~ /\.faa/}readdir(DIR); 
close DIR;

my %seqsH;
foreach (@files)
{
	open (IN, "<$anno_dir/$_"), or warn "cannot open $_\n";
	my @seqs = &gjoseqlib::read_fasta(\*IN);
	for my $i (0..$#seqs)
	{
		unless (($remove_x) && ($seqs[$i][2] =~ /x/i))
		{
			$seqsH{$seqs[$i][0]} = $seqs[$i][2]; 
		}
	}
	close IN;
}


my %hash; 
open (IN, "<$metadata"), or die "cannot open $metadata metadata file\n";
open (FAM, ">$prefix.fam"); 
open (NAM, ">$prefix.names"); 

while (<IN>)
{
	chomp; 
	s/Mature Envelope Glycoprotein Gn, Gn \(GP1\) Protein/Mature Envelope Glycoprotein Gn, Gn Protein/g; 
	s/Mature Envelope Glycoprotein Gc, Gc \(GP2\) Protein/Mature Envelope Glycoprotein Gc, Gc Protein/g; 
	my @array = split /\t/; 
	unless ($array[7] =~ /RNA/)
	{	
		my $anno = $array[-1]; 
		my $aa = $seqsH{$array[8]}; 
		my $name = $array[1]; 
		my $merge = $array[0];
		my $peg = $array[8]; 
		my $des = "$merge $name"; 
		push @{$hash{$anno}}, ([$peg, $des, $aa]);
		print FAM "$peg\t$array[2]\n"; 
		print NAM "$peg\t$name \[$merge\]\n"; 
	}
} 
close FAM;
close NAM;

my @fastas;
foreach (keys %hash)
{
	my $anno = $_; 
	$anno =~ s/.+\, //g; 
	$anno =~ s/ /_/g; 
	push @fastas, $anno;
	open (OUT, ">$prefix.$anno.aa"); 
	&gjoseqlib::print_alignment_as_fasta(\*OUT, @{$hash{$_}});
}

foreach (@fastas)
{
	system "mafft --thread 64 --reorder $prefix.$_.aa >$prefix.$_.fa";
	system "FastTree -wag <$prefix.$_.fa>$prefix.$_.nwk";
	system "remove_bad_tree_tips.pl <$prefix.$_.nwk >$prefix.$_.pruned.nwk";
	system "svr_tree_to_html -raw -a $prefix.names -c $prefix.fam <$prefix.$_.pruned.nwk>$prefix.$_.html";
}































