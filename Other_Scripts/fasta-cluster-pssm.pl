#! /usr/bin/env perl
use strict;
use gjoseqlib;
use Getopt::Long;
use Data::Dumper;
use Math::Round;
use Cwd;
use BlastInterface;


my $usage = 'fasta-cluster-pssm.pl [parms] <Protein-Sequences.fasta
			-h help
			-a "annotation string" (this also becomes the title in the pssm). 
			    If you select this program it will assume that you want all of your
			    sequences to be annotated the same. 
			-m min number of sequences for building pssm [d = 5]
			-f delete columns with less than frac bases (d = 0.33)
			-n delete columns from the n-terminus until >= frac identity (d = 0.75)
			-c delete columns from the c-terminus until >= frac identity (d = 0)
			-x include non-standard amino acids, by default it removes the 
			   entire sequence if it has ambiguous (B|J|X|Z) characters.
			-i include identical sequences [d = keep only unique]
			-d keep sequences beginning with a dash [d = throw away]
			-fd maximum fraction of dashes allowable for keeping a sequence (d = 0.20)
			-mi mmseqs identity [d = 0.8]
			-mc mmseqs coverage [d = 0.8]
			-k keep temp fasta (pre-aligned)
			-p "pssm prefix" (default = nothing).
			-tmp name temp fasta file 
';

my ($help, $tmp, $keep_x, $keep_i, $keep_tmp, $tmp, $start_dash, $annotation, $pssm_prefix);

my $min_seqs = 5;
my $f_insert = 0.33;
my $frac_n_term = 0.75;
my $frac_c_term = 0;
my $frac_dash = 0.2;
my $mmi = 0.8;
my $mmc = 0.8;


my $opts = GetOptions('h'     => \$help,
                      'a=s'   => \$annotation,
                      'm=s'   => \$min_seqs,
                      'f=f'   => \$f_insert,
                      'n=f'   => \$frac_n_term,
                      'c=f'   => \$frac_c_term,
                      'fd=f'  => \$frac_dash,
                      'x'     => \$keep_x,
                      'i'     => \$keep_i,
                      'k'     => \$keep_tmp,
                      'mmi=f' => \$mmi,
                      'mmc=f' => \$mmc,
                      'tmp=s' => \$tmp,
                      'p=s'   => \$pssm_prefix,
                      'd'     => \$start_dash);

if ($help){die "$usage\n"}
unless ($tmp){$tmp .= sprintf("%x", rand 16) for 1..10;}

## read proteins from stdin
my @prots = &gjoseqlib::read_fasta();

# create the tempdir, and go there.
my $base = getcwd;
mkdir ($tmp); 
chdir ($tmp);

### get rid of unique sequence, sequences with ambiguous residues, sort by length
my @prots2;
my $uniq = {};
if ($keep_i)
{	
	@prots2 = sort { length($b->[2]) <=> length($a->[2]) } @prots;
}
else
{
	if ($keep_x)
	{
		for my $i (0..$#prots)
		{
			$uniq->{$prots[$i][2]}->{ID}   = $prots[$i][0];
			$uniq->{$prots[$i][2]}->{ANNO} = $prots[$i][1];
		}
	}
	else
	{
		for my $i (0..$#prots)
		{
			unless ($prots[$i][2] =~ /(B|J|X|Z)/i)
			{	
				$uniq->{$prots[$i][2]}->{ID}   = $prots[$i][0];
				$uniq->{$prots[$i][2]}->{ANNO} = $prots[$i][1];
			}
		}
	}

	my @unsorted;
	foreach (keys %$uniq)
	{
		my $seq  = $_;
		my $anno = $uniq->{$_}->{ANNO};
		my $id   = $uniq->{$_}->{ID}; 
		push @unsorted, [$id, $anno, $seq];
		@prots2 = sort { length($b->[2]) <=> length($a->[2]) } @unsorted;
	}
}

# Print fasta to temp
open (OUT, ">$tmp.fasta");

for my $i (0..$#prots2)
{
	my $id = $prots2[$i][0];
	$id =~ s/\|[^|]*$// if $id =~ tr/|// > 1;
	if ($annotation)
	{
		&gjoseqlib::print_alignment_as_fasta(\*OUT, ([$id, $annotation, $prots2[$i][2]]));
	}
	else
	{	
		&gjoseqlib::print_alignment_as_fasta(\*OUT, ([$id, $prots2[$i][1], $prots2[$i][2]]));
	}
}


print STDERR "Running mmseqs\n"; 
# Run mmseqs
system "mmseqs easy-cluster $tmp.fasta $tmp.mmseq tmp --min-seq-id $mmi -c $mmc --cov-mode 0 >/dev/null";




#Process the clusters into individual fasta files.
 
open (IN, "<$tmp.mmseq_all_seqs.fasta"), or die "cannot open all_seqs fasta file from mmseqs\n"; 
my $seqH = {};
my @seqs = &gjoseqlib::read_fasta(\*IN); # garys read program ignores fasta identifiers with out sequence.
close IN;
for my $i (0..$#seqs)
{
	$seqH->{$seqs[$i][0]}->{ANNO} = $seqs[$i][1];
	$seqH->{$seqs[$i][0]}->{SEQ}  = $seqs[$i][2];
}

my %clusters;
open (IN, "<$tmp.mmseq_cluster.tsv"), or die "couldn't open tsv file from mmseqs\n"; 
while (<IN>)
{
	chomp;
	my ($id, $memb) = split /\t/;
	push @{$clusters{$id}}, $memb;
}
close IN;

my $count = 1;
mkdir ("clusters"); 
open (OUT2, ">Leftover_Seqs.aa");

foreach (sort { scalar(@{$clusters{$b}}) <=> scalar(@{$clusters{$a}}) } keys %clusters)
{	
	my @array = @{$clusters{$_}};
	
	if ((scalar @array) > $min_seqs)
	{
		open (OUT, ">clusters/$count.fasta");
		foreach (@array)
		{
			&gjoseqlib::print_alignment_as_fasta(\*OUT, ([$_, $seqH->{$_}->{ANNO}, $seqH->{$_}->{SEQ}]));
		}
		close OUT;
		$count ++;
	}
	else
	{
		foreach (@array)
		{
			&gjoseqlib::print_alignment_as_fasta(\*OUT2, ([$_, $seqH->{$_}->{ANNO}, $seqH->{$_}->{SEQ}]));
		}
	}
}
close OUT2;
system "mafft --thread 24 --quiet --reorder $base/$tmp/Leftover_Seqs.aa>Leftover_Seqs.fa";



mkdir ("alis"); 
mkdir ("corrected_alis"); 
mkdir ("pssms"); 

opendir (DIR, "$base/$tmp/clusters");
my @clusterF = grep{$_ !~ /^\./}readdir(DIR); 
close DIR;

open (REP, ">Curation_Report");
print REP "Alignment\tAli_Len\tFirst_AA\tOriginal Seqs\tFinal Seqs\tLow_Occ_Cols_Cut\tN_Term_Cut\tC_Term_Cut\tDash_Starts_Removed\tCut_Cols\n";

foreach (@clusterF)
{
	my $file = $_;
	my $ali_file = $_;
	$ali_file =~ s/\.fasta//g;
	print STDERR "\n\nAligning $file\n"; 
	open (IN, "mafft --thread 24 --quiet --reorder $base/$tmp/clusters/$file |");
	my @ali = &gjoseqlib::read_fasta(\*IN); 
	close IN; 

	#print the original alignment
	#I might eventually want to turn this off.
	open (OUT, ">alis/$ali_file.fa"); 
	&gjoseqlib::print_alignment_as_fasta(\*OUT, @ali);
	close OUT;
	
	
	# get rid of any sequences containing too many dashes.
	# 
	my @ali2;
	for my $i (0..$#ali)
	{
		my $len = length ($ali[$i][2]);
		my $dashes = $ali[$i][2] =~ tr/-//;

		if (($dashes/$len) <= $frac_dash) 
		{
			push @ali2, ($ali[$i]);
		}
	}

	my $n_seqs_orig = scalar @ali;
	my $n_seqs_internal_dash = scalar @ali2;
	
#	print STDERR "Dash removal Before = $n_seqs_orig\tAfter = $n_seqs_internal_dash\n"; 
	
	# Now I need to pack the alignment, this removes empty columns.
	my @ali3 = &gjoseqlib::pack_alignment(@ali2);
	

	
	## Identify alignment inserts < $f_insert

	my $n_seqs = scalar @ali3;
	my $min_n = (round ($f_insert * $n_seqs)); 
#	print STDERR "$file\tTotal_SEQs = $n_seqs\tMin Seqs = $min_n\n"; 
	my %col_sum;     #number non dash characters
	my $aa_sum = {};   #aa count
	
	for my $i (0..$#ali3)  #count the number of characters in the column.
	{
		my $sequence = uc $ali3[$i][2];
		my @bases = split ("", $sequence);
		for my $j (0..$#bases)
		{
			if ($bases[$j] =~ /\w/)
			{
				$aa_sum->{$j}->{$bases[$j]} ++;  #only counts letters not dashes
				$col_sum{$j} ++;
			}
		}
	}

	## Add low occupancy columns to the exclusion hash 
	my %exclude;
	foreach (keys %col_sum)
	{
		my $f_occ = ($col_sum{$_}/$n_seqs);
		if ( $f_occ < $f_insert)
		{
			print STDERR "Excluding $_\t$f_occ\n"; 
			$exclude{$_} = 1;
		}
	}	
	my $low_occ = keys %exclude; #count of the low occupancy columns removed.
	my $ali_len = length ($ali3[0][2]);
	
	# Trim from the N-terminus
	my $n_term_trimmed = 0;
	my $c_term_trimmed = 0;

	for my $k (0..$ali_len)
	{
		my $aa_resR = $aa_sum->{$k};
		my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys (%$aa_resR);
		my $mc_res = $sorted_res[0];
		my $cons = ($aa_sum->{$k}->{$mc_res}/$n_seqs_internal_dash);

		if ($cons < $frac_n_term)
		{
			$exclude{$k} = 1;
			$n_term_trimmed ++;
			print STDERR "Trimmed $k $mc_res, $aa_sum->{$k}->{$mc_res}, $cons, from N-terminal\n"
		}
		else
		{
			last;
		}
	}
	
	# trim backwards from C-terminus
	for (my $l = ($ali_len -1); $l >= 0; $l --)
	{
		my $aa_resR = $aa_sum->{$l};
		my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys (%$aa_resR);
		my $mc_res = $sorted_res[0];
		my $cons = ($aa_sum->{$l}->{$mc_res}/$n_seqs_internal_dash);

		if ($cons < $frac_c_term)
		{
			$exclude{$l} = 1;
			$c_term_trimmed ++;
			print STDERR "Trimmed $l, $mc_res, $aa_sum->{$l}->{$mc_res},$cons, from C-terminal\n"
		}
		else
		{
			last;
		}	
	}
	
	
	# read through the alignment, unless the position exists in $exclude, add it
	# Get rid of sequences starting with a dash, unless -d is turned on
	my @ali4;
	my $dash_starts_removed = 0;
	for my $i (0..$#ali3)
	{
		my $id = $ali3[$i][0];
		my $anno = $ali3[$i][1];
		my @bases = split ("", $ali3[$i][2]);
		my $string;
		for my $j (0..$#bases)
		{	
			unless (exists $exclude{$j})
			{
				$string .= $bases[$j]; 
			}	
		}
		
		if ($start_dash)
		{
			push @ali4, ([$id, $anno, $string]); 
		}
		elsif ((! $start_dash) && ($string !~ /^-/))
		{
			push @ali4, ([$id, $anno, $string]); 
		}
		elsif ((! $start_dash) && ($string =~ /^-/))
		{
			$dash_starts_removed ++;
		}
	
	}
	
	my @ali5 = &gjoseqlib::pack_alignment(@ali4);
	
	
	# Get rid of identical sequences
	my $unique = {};
	for my $i (0..$#ali5)
	{
		my $seq = uc $ali5[$i][2];
		$unique->{$seq}->{ID}   = $ali5[$i][0];
		$unique->{$seq}->{ANNO} = $ali5[$i][1];
	}
	my @ali6;
	#my @ali7;  # I need internally hide the vertical bar to make a psiblast-compatible id.
	foreach (keys %$unique)
	{
		push @ali6, ([$unique->{$_}->{ID}, $unique->{$_}->{ANNO}, $_]); 
		my $id = $unique->{$_}->{ID};
		
		#$id =~ s/\|/##/g; 
		#push @ali7, ([$id, $unique->{$_}->{ANNO}, $_]); 

	}

	my $n_seqs_final = scalar @ali6;
	my $start_char =  substr($ali6[0][2], 0, 1);
	my $ali_len = length $ali6[0][2];

	print STDERR  "Original_Seqs = $n_seqs_orig,  Final_Seqs = $n_seqs_final\n";
	print STDERR "Low Occ Cols Cut = $low_occ, N-term Cut = $n_term_trimmed, C-term Cut = $c_term_trimmed\tDash Starts Removed = $dash_starts_removed\n\n";

	print REP "$ali_file\t$ali_len\t$start_char\t$n_seqs_orig\t$n_seqs_final\t$low_occ\t$n_term_trimmed\t$c_term_trimmed\t$dash_starts_removed\t";
	my @sorted_exclude = sort { $a <=> $b } keys %exclude;
	print REP join (",", @sorted_exclude), "\n"; 
	
	open (OUT, ">corrected_alis/$ali_file.fa"); 
	&gjoseqlib::print_alignment_as_fasta(\*OUT, @ali6);
	close OUT;
	
	if ($n_seqs_final >= $min_seqs)
	{

		my $pssm_file = $ali_file;
		if ($pssm_prefix){$pssm_file = "$pssm_prefix.$ali_file";}
	
		my %pssm_opts = (outPSSM  => "pssms/$tmp.pssm");
		my $pssm = BlastInterface::alignment_to_pssm(\@ali6, \%pssm_opts);
	
		open (IN, "<pssms/$tmp.pssm");
		open (OUT, ">pssms/$pssm_file.pssm");
		
		my $des = 0;

		#blast interface is useful for making the pssm, but it doesn't format the title
		#nicely, so below i loop over that, and fixed it.  Basically, it enumerates 
		#each title, which is ugly, and if you don't give it a title, it will give you an 
		#undefined title plus your annotation string, so I just pick it up from the 
		#annotation string, and delete the undefined one that is printed compulsively.
		
		while (<IN>)
		{
			chomp;
			if (($des > 1) && ($_ =~ /title/))
			{
				next;
			}
			elsif (($des > 1) && ($_ =~ /\}/))
			{
				$des = 0;
				next;
			}
			elsif ($_ =~ /descr \{/)
			{
				$des ++;
				if ($des <=1)
				{
					print OUT "$_\n";
				}
			}
			elsif ($des <= 1)
			{
				print OUT "$_\n";
			}
		}		
		unlink ("pssms/$tmp.pssm");	
	}
}

close REP;




























