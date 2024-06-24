#! /usr/bin/env perl
use strict;
use gjonewicklib;
use Data::Dumper;
use Getopt::Long;
use gjostat;
use ffxtree;


my $usage = 'remove_bad_tree_tips.pl -n (n * stdev) <tree.nwk > clean_tree.nwk

	This program reads in a tree file, collapses identical sequences, computes
	an average tip distance, and then removes long branches that are (N * STDEV) 
	longer than the average tip distance.  It returns the tree with out the long tips.
	
	-n (multiplier to the stdev d = 1.25)
	-h help



';

my ($n, $help,);
my $opts = GetOptions('h'      => \$help,
                      'n=f'    => \$n);


if ($help){die "$usage\n";}
unless ($n){$n = 1.25}; 

my $TREE = \*STDIN; 
my $treestr = join( "", <$TREE> ); 
my $tree1 = &gjonewicklib::parse_newick_tree_str( $treestr ); 
my @tips = newick_tip_list ($tree1);

#doesn't work in gjonewicklib;
#my $tree2 = &gjonewicklib::collapse_zero_length_branches( $tree1 );

# compute an upper bound for the distance based on an average an stdev
# computed off of all the zero-length branches.  

my $tree2 = &ffxtree::collapse_identical_seqs($tree1);
my %tip_dists2 = newick_tip_distances( $tree2 );
my @dists = values %tip_dists2;
my ($mean, $stdev) = &gjostat::mean_stddev(@dists);
my $upper = ($mean + ($stdev * $n));
#my $lower = ($mean - ($stdev * 1.25));  #lower doesn't work right


# purge the original ids based on the upper bound of the distance
my %tip_dists1 = newick_tip_distances( $tree1 );
my %exclude;
my @keep_tips;

foreach (sort {$tip_dists1{$a} <=> $tip_dists1{$b}} keys %tip_dists1)
{
	if ($tip_dists1{$_} > $upper)
	{
		print STDERR "$_\t$tip_dists1{$_}\tSTDEV=$stdev\tUPPER=$upper\n";
	}
	else
	{
		push @keep_tips, $_;
	}
}
#my $newtree = &gjonewicklib::newick_subtree( $tree1, @keep_tips); 
my $newtree = &gjonewicklib::newick_subtree( $tree1, @keep_tips); 
&gjonewicklib::writeNewickTree($newtree); 




# this didn't work out.
#95% CI
#my $n = scalar @dists; 
#my $se = $stdev/sqrt($n); 
#my $z = 2.576; 
#my $me = ($z  * $se); 
#my $lower = ($mean - $me);
#my $upper = ($mean + $me);





















