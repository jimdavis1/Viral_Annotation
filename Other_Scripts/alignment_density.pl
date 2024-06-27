#! /usr/bin/env perl
use strict; 
use Getopt::Long;
use gjoseqlib;
use gjostat;
use Math::Round;
use Data::Dumper;

my $usage = 'cut_ali_inserts.pl [parms] <Sequence.fasta
			
			This program computes a density score for each sequence in the alignment
			that is defined in the following way:
			
			1.  The fraction occupancy of each column is computed (#characters/#dashes)
			2.  The occupancy scores for each  occupied column are summed 
			    for each sequence.
			3.  Sequences that are >= -z (d = 4) stdevs from the mean are removed.
				
			4. it will also remove any sequence whose average occupancy is >= z away 
			   from the average of all sequences	
			
			5. it will also remove any sequence that is >= -z stdevs from the average length.
						
			-h help
			-p prefix for outputfile names
			-z z value cutoff D = 4
			-t print the table of scores (ID, Z, Density, Avg Des, Stdev, DENS)
			-m min protein length (d = 100)

';
my $min_len = 100;
my $cutoff = 4;
my ($help, $prefix, $table);
my $opts = GetOptions('h'   => \$help,
                      'p=s' => \$prefix,
                      't'   => \$table,
                      'm=i' => \$min_len,
                      'z=f' => \$cutoff);

if ($help){die "$usage\n"}
unless ($prefix){die "must declare a -p prefix for output files\n"};


my @seqs = &gjoseqlib::read_fasta(); 


my %col_sum;
my $total_seqs = scalar @seqs;
my %length;

for my $i (0..$#seqs)
{
	my @bases = split ("", $seqs[$i][2]);

	for my $j (0..$#bases)
	{
		if ($bases[$j] =~ /\w/)
		{
			$col_sum{$j} ++;
			$length{$seqs[$i][0]} ++;
		}
	}
}

my %frac_occ;
foreach (keys %col_sum)
{
	$frac_occ{$_} = ($col_sum{$_}/$total_seqs);
}



my @all_mean;
my @all_dens;
my $seq_data = {};

# compute the sequence density and avg density for each sequence
for my $i (0..$#seqs)
{
	my @bases = split ("", $seqs[$i][2]);
	my @array;
	my $sum;
	for my $j (0..$#bases)
	{
		if ($bases[$j] =~ /\w/)
		{
			push @array, $frac_occ{$j};
			$sum += $frac_occ{$j};
		}
	}
	my ($mean, $stdev) = &gjostat::mean_stddev(@array);
	
	$seq_data->{$seqs[$i][0]}->{DES}  = $seqs[$i][1];
	$seq_data->{$seqs[$i][0]}->{SEQ}  = $seqs[$i][2];
	$seq_data->{$seqs[$i][0]}->{DENS} = $sum;
	$seq_data->{$seqs[$i][0]}->{AVG}  = $mean;

	push @all_mean, $mean;
	push @all_dens, $sum;
	
	#print "$seqs[$i][0]\t$sum\t$mean\t$stdev\t$seqs[$i][1]\n";
}

#compute an avg for density sums across sequences
my ($mean_dens, $stdev_dens) = &gjostat::mean_stddev(@all_dens);

#compute position avg per sequence
my ($mean_pos, $stdev_pos) = &gjostat::mean_stddev(@all_mean);

#compute avg length per sequence
my ($mean_len, $stdev_len) = &gjostat::mean_stddev(values %length);



open (YES, ">$prefix.keep.fa");
open (NO, ">$prefix.cull.fa");

my @keep;
foreach (keys %$seq_data)
{
	my $id = $_;
	my $z_dens = abs(($seq_data->{$id}->{DENS} - $mean_dens)/$stdev_dens);
	my $z_pos  = abs(($seq_data->{$id}->{AVG} - $mean_pos)/$stdev_pos);
	my $z_len  = abs(($length{$id} - $mean_len)/$stdev_len);

	if ($table)
	{
		#print  "$id\t$\t$seq_data->{$id}->{DENS}\t$mean_dens\t$stdev_dens\n";
		print "$id\t$z_dens\t$z_pos\t$z_len\t$seq_data->{$id}->{DENS}\t$seq_data->{$id}->{AVG}\t$length{$id}\n";
	
	}
	if (($z_dens < $cutoff) && ($z_pos < $cutoff) && ($z_len < $cutoff) && ($length{$id} >= $min_len))
	{
		push @keep, ([$id, $seq_data->{$id}->{DES}, $seq_data->{$id}->{SEQ}]); 
	}
	else
	{
		&gjoseqlib::print_alignment_as_fasta(\*NO, [$id, $seq_data->{$id}->{DES}, $seq_data->{$id}->{SEQ}]); 
	}	
}

my @packed = pack_alignment(@keep);
&gjoseqlib::print_alignment_as_fasta(\*YES, @packed);
close YES;
close NO;












