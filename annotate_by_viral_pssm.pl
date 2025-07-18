#! /usr/bin/env perl

use strict;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;
use Cwd;
use gjoseqlib;

my $default_data_dir = $ENV{LOWVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

my $usage = 'annotate_by_viral_pssm.pl [options] -i subject_contig(s).fasta 

		-h   help
		-i   Input subject contigs in fasta format
		-t   Declare a temp file (d = random)
		-tax Declare a taxonomy id (D = 10239)
		-g   Genome name (D = Viruses);
		-p   Prefix for the output files (.ffn, .faa, and .tbl files)
		-s   Append sequences to the feature table (default is that they are left off)
		-ks  If a pssm match extends beyond a stop codon it will keep the entire region of the match including the stop. 
		-threads number of blast threads in the blastn and tblastn

		-min Minimum contig	length (d = 1000)  # otherwise the genome is rejected
		-max Maximum contig length (d = 30000) # for reference Measles is 15,894 and Beilong is 19,212
	
		-mcb Minimum contig bitscore to enable annotation (d = 500) #otherwise the genome is rejected.

        -j   Full path to the options file in JSON format which carries data for a match (D = $default_data_dir/Viral_PSSM.json)
		-c   Representative contigs directory (D = $default_data_dir/Viral-Rep-Contigs)
		-pssm   Base directory of PSSMs   (D = $default_data_dir/Viral-PSSMs)
	           Note that this is set up as a directory of pssms
	           right now this is hardcoded as: "virus".pssms within this directory.

	      
	   Debug parms (turns off output file generation)
	   	-tmp keep temp dir
		-no  no output files generated, to be used in conjunction with one of the following:
		    -dna print only genes to STDOUT 
		    -aa print proteins to STDOUT
		    -tbl print only feature table to STDOUT
		    -ctbl [file name] concatenate table results to a file (for use with many genomes)


 	## Output Table columns are the following:
 			0	genome_id	
 			1	genome_name	
 			2	accession (original contig id)
 			3	annotation_source (always LV)	
 			4	feature_type	
 			5	a made up feature id	
 			6	gene symbol	
 			7	start	
 			8	end	
 			9	strand	
 			10	na_length	
 			11	annotation taxon
 			12	closest genome
 			13	closest genome bit score (BLASTn of contigs)
 			14  closest genome id(s)
 			15  closest genome name
 			16	PSSM
 			17	NT sequence (if -s append sequences)	
 			18	AA sequence (if -s append sequences)
';

my ($help, $opt_file, $contig_file, $tmp, $tax, $keep_stop, $genome_name, $cdir, $pdir, $keep_temp, $min_len, $max_len, $aa_only, $dna_only, $tbl_only, $no_out, $ctbl, $prefix, $append_seqs, $threads, $min_contig_bit);

my $opts = GetOptions( 'h'         => \$help,
                       'tmp'       => \$keep_temp,
                       'no'        => \$no_out,
                       'aa'        => \$aa_only,
                       'dna'       => \$dna_only,
                       'tbl'       => \$tbl_only,
                       'ctbl=s'    => \$ctbl,
                       'i=s'       => \$contig_file,
                       'tax=s'     => \$tax,
                       'threads=s' => \$threads,
                       't=s'       => \$tmp,
                       'g=s'       => \$genome_name,
                       'c=s'       => \$cdir,
                       'pssm=s'    => \$pdir,
                       'min=s'     => \$min_len,
                       'max=s'     => \$max_len,
                       'mcb=f'     => \$min_contig_bit,
                       'j=s'       => \$opt_file,
                       'ks'        => \$keep_stop,
                       'p=s'       => \$prefix,
                       's'         => \$append_seqs); 

if ($help){die "$usage\n";}
unless ($contig_file ){die "must declare an input subject file with -i \n\n$usage\n";}

# Name the Temp file:
# generates a random 20 digit string of 0-9a-f
unless ($tmp){$tmp .= sprintf("%x", rand 16) for 1..20;}
unless ($min_len){$min_len = 1000; }
unless ($max_len){$max_len = 30000; }
unless ($min_contig_bit){$min_contig_bit = 500;}
unless ($tax){$tax = "10239"; }
unless ($genome_name){$genome_name = "Viruses"; }
unless ($prefix){$prefix = "Viral_Annotation";}
unless ($cdir){$cdir = "$default_data_dir/Viral-Rep-Contigs"; }
unless ($pdir){$pdir = "$default_data_dir/Viral-PSSMs"; }
unless ($opt_file){$opt_file = "$default_data_dir/Viral_PSSM.json"; }
unless ($threads){$threads = 8}


## This is the json file with all of the protein-specific information for 
## Customizing the annotation.
open (IN, "<$opt_file"), or die "Cant find JSON options file, use -opt\n"; 
my $options = decode_json(scalar read_file(\*IN));
close IN;


# Version the Genome (another random alpha numeric)
my $version;
{$version .= sprintf("%x", rand 16) for 1..6;}  ## 11M combos

# Open file handles for the outputs

unless($no_out)
{
	if ($prefix)
	{
		open (AA,  ">$prefix.faa");
		open (DNA, ">$prefix.ffn");
		open (TBL, ">$prefix.feature.tbl");
	}
	else
	{
		open (AA, ">$tax.$version.faa");
		open (DNA, ">$tax.$version.ffn");
		open (TBL, ">$tax.$version.feature.tbl");
	}
}

if ($ctbl)
{
	open(CTBL, ">>", "$ctbl")  
}

# Make a hash out of the subject contigs
open (IN, "<$contig_file");
my @seqs = &gjoseqlib::read_fasta(\*IN);
close IN;

# Make sure it falls within the min and max contig length.
# Count up the size so that it can handle multiple contigs.
my $len = 0;
my %contigH;
my @contig_order;
for my $i (0..$#seqs)
{
	$contigH{$seqs[$i][0]} = uc ($seqs[$i][2]);
	$len += length($seqs[$i][2]);
	push @contig_order, $seqs[$i][0];
}
if (($len < $min_len) || ($len > $max_len))
{
	die "Unexpected input genome size equal to: $len\nMin = $min_len\tMax = $max_len\n";	
}	

# Make the temp dir.
my $base = getcwd;
mkdir ($tmp); 
system "cp $contig_file $tmp";
my $s_file = $contig_file;
$s_file =~ s/.+\///g; 

chdir ($tmp);
system "makeblastdb -dbtype nucl -in $s_file >/dev/null";


#   Select the correct species based on the contig blastN
#   If this begins to break down, new reference contigs can be added to the directory
#   Formatting them as genus.version.dna

opendir (DIR, "$cdir");
my @reps = grep{$_ !~ /^\./}readdir(DIR); 
close DIR;

my $best_contig_bit = 0;
my $best_virus_match;
my $best_rep;
foreach (@reps)
{
	my $rep = $_;
	my $virus = $_;
	$virus =~ s/\..+//g; 
	
	my $rep_file = "$cdir/$rep";
	open (IN, "blastn -query $rep_file -subject $s_file -evalue 0.5 -reward 2 -penalty -3 -word_size 11 -outfmt 15 -soft_masking false -num_threads $threads 2>/dev/null |") or die "Could not run blastn: $!"; 
	my $blastn = decode_json(scalar read_file(\*IN));	
	close IN;

	my $match_bit = get_blastn_bit($blastn);  
	unless ($match_bit < $best_contig_bit)
	{
		$best_contig_bit = $match_bit;
		$best_virus_match = $virus;
		$best_rep = $rep;
	}
}	

# here we make it fail gracefully if we don't enounter a reference genome with a high 
# enough blastn bit score.

if ($best_contig_bit < $min_contig_bit) {
    print STDERR "No matching reference contigs with bit score greater than $min_contig_bit\n";
    exit(1);
}



#get closest genome data
my $best_rep_ids   = $options->{$best_virus_match}->{close_genomes}->{$best_rep}->{genome_ids};
my $best_rep_name  = $options->{$best_virus_match}->{close_genomes}->{$best_rep}->{genome_name};







my $virus = $best_virus_match;	
print STDERR "-----------------------\nMatching $virus\tBit = $best_contig_bit\n-----------------------\n"; 
opendir (DIR, "$pdir/$virus.pssms");
my @pssm_dirs = grep{$_ !~ /^\./}readdir(DIR);  # reads the directory of PSSM dirs
close DIR;

my @all_seqs;
my %positions;
#my ($upstream_ext, $downstream_ext, $keep_stop, $bit_cutoff, $cov_cutoff);

my $non_pssm_feat = {};
foreach (@pssm_dirs)  #Each PSSM dir contains one or more PSSMs for a given homolog
{
	my $pssmdir = $_; 
	opendir (DIR, "$pdir/$virus.pssms/$pssmdir");
	my @pssm_files = grep{$_ !~ /^\./}readdir(DIR); # gets the pssms
	close DIR;

	# retrieve the protein specific options for the corresponding PSSM.
	unless ($options->{$virus}->{features}->{$pssmdir})
	{
		print STDERR "No data in JSON file for VIRUS: $virus\tPSSM: $pssmdir\n"; 
		next;
	}
	my $upstream_ext   = $options->{$virus}->{features}->{$pssmdir}->{upstream_ext};
	my $downstream_ext = $options->{$virus}->{features}->{$pssmdir}->{downstream_ext};
	my $bit_cutoff     = $options->{$virus}->{features}->{$pssmdir}->{bit_cutoff};
	my $cov_cutoff     = $options->{$virus}->{features}->{$pssmdir}->{coverage_cutoff};
	my $start_to_met   = $options->{$virus}->{features}->{$pssmdir}->{start_to_met};
	my $feature_type   = $options->{$virus}->{features}->{$pssmdir}->{feature_type};
	my $anno           = $options->{$virus}->{features}->{$pssmdir}->{anno};
	my $symbol         = $options->{$virus}->{features}->{$pssmdir}->{gene_symbol};
		
	print STDERR "\t$virus\t$pssmdir\t$anno\tbit\t$bit_cutoff\tcov\t$cov_cutoff\tkeep_stop\t$keep_stop\tupstream_ext\t$upstream_ext\tdownstream_ext\t$downstream_ext\n"; 		

	#   Select the best pssm per protein
	#   If this begins to break down, new reference pssms can be added to the 
	#   "Taxon.pssms/protein/"   directory
	#   Adding more pssms should theoretically improve the accuracy

	my ($best_pssm, $best_results, $pssm_name);
	my $best_bit = 0;

	foreach (@pssm_files)  ##  <---- This is where I would need to lookup the rules for each protein
	{		
		my $pssm_file = "$pdir/$virus.pssms/$pssmdir/$_";
		my $name = $_;
		open (IN, "tblastn -outfmt 15 -db $s_file -in_pssm $pssm_file -num_threads $threads |") or die "Could not run tblastn: $!";
		my $pssm_blast = decode_json(scalar read_file(\*IN));		
		close IN; 	
		
		my  ($results, $hsp_best_bit) = matching_tblastn_hsps_json($pssm_blast, $bit_cutoff, $cov_cutoff);
		print STDERR "\t$name\t$hsp_best_bit\n"; 

		unless ($hsp_best_bit < $best_bit)   #evaulate each pssm blast based on the bit score, and pick the best one
		{
			$best_bit = $hsp_best_bit;
			$best_pssm = $name;
			$best_results = $results;
		}
	}
	print STDERR "\n\tChoosing $best_pssm\t$best_bit\n\n"; 
	my $nhsps = scalar @$best_results;
	for my $i (0..($nhsps -1))
	{		
		my $contig = $best_results->[$i]->{contig};
		$contig =~ s/\s.+//g;								
		my $hseq = $best_results->[$i]->{hseq};
		$hseq =~ s/\-//g; # eliminate alignment gaps in hit seq
			
		print STDERR "\t$contig\t$anno\t$best_results->[$i]->{hit_from}\t$best_results->[$i]->{hit_to}\t$best_results->[$i]->{frame}\t$best_results->[$i]->{bit}\n";

		## Okay, Now we need to find our features.
	
	
		# lookup DNA Coordinates
		my ($from, $to);
		my $strand = "+";
		if ($best_results->[$i]->{frame} < 0)
		{
			$from = $best_results->[$i]->{hit_to};
			$to   = $best_results->[$i]->{hit_from};
			$strand = "-";
		}
		else
		{
			$from = $best_results->[$i]->{hit_from};
			$to   = $best_results->[$i]->{hit_to};
		}

		###  $from and $to are the coordinates of the protein from the blast without the stop codon.
		###  The end extension loop below, scans in the 3' direction for the next stop codon.
		###  if $downstream_ext is declared,  it will get new gene coordinates with the stop
		###  codon being included in the gene coordinates by convention. 		

		my ($gene_begin, $gene_end);
		if ((! $downstream_ext) || ($hseq =~ /\*/))
		{
			$gene_begin = $from;
			$gene_end = $to;
		}
		else
		{
			($gene_begin, $gene_end)= scan_to_stop_codon ($from, $to, $contigH{$contig});
			print STDERR "\tScanning to Stop codon, Original Coords: $from to $to\tNew Coords: $gene_begin to $gene_end\n"; 	

		}
		#print STDERR "COV = $best_results->[$i]->{cov}\n";

		
		## If a match contains an internal stop, everything after the stop is cropped by default.		
		my $new_end;
		if (($hseq =~ /\*/) && (! $keep_stop))
		{
			print STDERR "\tMatch contains a stop codon, cropping\n";
			$new_end = crop_to_stop_codon($gene_begin, $gene_end, $hseq);
		
			#check the coverage.
			my $cov2 =($len/(abs(($best_results->[$i]->{q_to} - $best_results->[$i]->{q_from})+1))); 

			if ($cov2 < $cov_cutoff)
			{
				print STDERR "\tCoverage prior to stop codon \($cov2\) is lower than cutoff: $cov_cutoff\n";
				next;
			}
			print STDERR "\tCropping to coordinates: $gene_begin\tto\t$new_end\n"; 			
		}	
		if ($new_end){$gene_end = $new_end};	
		

		# When the pssm match doesn't start with an M, it will search to the left until it finds the first one.
		# This must be turned off if there is a non-AUG start.

		my ($new_gene_begin, $new_gene_end); 
		if (($hseq !~ /^M/i) && ($upstream_ext))
		{
			($new_gene_begin, $new_gene_end) = scan_to_met_start( $gene_begin, $gene_end, $contigH{$contig});
			print STDERR "\tScanning for Met start, Original Coords: $gene_begin to $gene_end\tNew Coords: $new_gene_begin to $new_gene_end\n";

		}
		if ($new_gene_begin){$gene_begin = $new_gene_begin};	


		my $gene = &gjoseqlib::DNA_subseq($contigH{$contig}, $gene_begin, $gene_end );	
		my $protein = &gjoseqlib::translate_seq( $gene );

		if ($start_to_met){$protein =~ s/^[A-Z]/m/i;}
	
		# Set up calling non-pssm features that are anchored to pssm coordinates
		if (exists $options->{$virus}->{features}->{$pssmdir}->{non_pssm_partner})
		{
			my @non_pssms = @{$options->{$virus}->{features}->{$pssmdir}->{non_pssm_partner}};
			
			for my $i (0..$#non_pssms)
			{
				my $start_coord;
				my $stop_coord;
				my $feat = $non_pssms[$i];
				#if the current pssm match feature corresponds with the non-pssm start site
				if ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{START_OFFSET} = $options->{$virus}->{features}->{$feat}->{begin}->{begin_offset};
					if ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm_loc} =~ /START/)
					{	
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_begin;
					}
					elsif ($options->{$virus}->{features}->{$feat}->{begin}->{begin_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_end;
					}
				}
				elsif ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{STOP_OFFSET} = $options->{$virus}->{features}->{$feat}->{end}->{end_offset};
					if ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm_loc} =~ /START/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_begin;
					}
					elsif ($options->{$virus}->{features}->{$feat}->{end}->{end_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_end;
					}
				}			
				$non_pssm_feat->{$feat}->{ANNO}     = $options->{$virus}->{features}->{$feat}->{anno};
				$non_pssm_feat->{$feat}->{MIN}      = $options->{$virus}->{features}->{$feat}->{min_len};
				$non_pssm_feat->{$feat}->{MAX}      = $options->{$virus}->{features}->{$feat}->{max_len};			
				$non_pssm_feat->{$feat}->{AA}       = $options->{$virus}->{features}->{$feat}->{translate};
				$non_pssm_feat->{$feat}->{TYPE}     = $options->{$virus}->{features}->{$feat}->{feature_type};
				$non_pssm_feat->{$feat}->{SYMBOL}   = $options->{$virus}->{features}->{$feat}->{gene_symbol};
			}
		}
			push @all_seqs, ([$best_results->[$i]->{contig}, $gene_begin, $gene_end, $anno, $strand, $best_pssm, $gene, $protein, $feature_type, $symbol]); 
	}
	print STDERR "-----------------------\n"; 
}


#add any non-pssm features that are anchored to PSSM coordinates
if ($non_pssm_feat)
{
	my @tuples = call_non_pssm_features($non_pssm_feat, %contigH);
	push @all_seqs, @tuples;
}


#Tuple is:	
# 0 contig
# 1 start
# 2 end
# 3 anno
# 4 strand
# 5 pssm  (this is the full pssmfile name)
# 6 gene
# 7 protein
# 8 feature_type	
# 9 symbol

# $virus is the best virus match
# Sort the output in order of contig and then start position. 

my (@prot_seqs, @gene_seqs);
my $count = 0;
foreach (@contig_order)
{
	my $contig = $_;
	my @features = ();
	for my $i (0..$#all_seqs)
	{
		if ($all_seqs[$i][0] eq $contig)
		{
			push @features, $all_seqs[$i];		
		}		
	}
	my @sorted =  sort { $a->[1] <=> $b->[1] } @features;
	foreach(@sorted)
	{
		$count ++; 
		# this is just a silly place holder unitl bob can help me get set up with official
		# ID generation
		my $prot_id = "lv\|$tax\.$version\.$_->[8]\.$count"; 
		
		#Trimming the star from the end of the protein sequence. 
		#This could be removed if its unwanted.
		my $protein = $_->[7];
		$protein =~ s/\*$//g;

		push @gene_seqs, ([$prot_id, $_->[3], $_->[6]]);
		push @prot_seqs, ([$prot_id, $_->[3], $protein]);
		
		# Original Feature Table looks like this:
		# genome_id	genome_name	accession	annotation_source	feature_type	patric_id	refseq_locus_tag	start	end	strand	na_length	gene	product	plfam_id	pgfam_id
		# For now, I will keep:
		# genome_id	genome_name	accession(contig)	annotation_source	feature_type	madeupid	Gene_symbol	start	end	strand	na_length	Fam	Closest_Genome	best_contig_bit pssm	product	NT_Seq	AA_seq


		my $na_len = length $_->[6];
		unless($no_out)
		{
			if ($append_seqs)
			{
				print TBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";
			}
			else
			{
				print TBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";
			}
		}
		if($tbl_only)
		{
			if ($append_seqs)
			{
				print "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";
			}
			else
			{
				print "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";
			}		
		}
		if($ctbl)
		{
			if ($append_seqs)
			{
				print CTBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\t$_->[6]\t$_->[7]\n";

			}
			else
			{
				print CTBL "$tax\.$version\t$genome_name\t$_->[0]\tLV\t$_->[8]\t$prot_id\t$_->[9]\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$virus\t$best_rep\t$best_contig_bit\t$best_rep_ids\t$best_rep_name\t$_->[5]\t$_->[3]\n";

			}
		}
	}
}




unless($no_out)
{
	&gjoseqlib::print_alignment_as_fasta(\*DNA, @gene_seqs);
	&gjoseqlib::print_alignment_as_fasta(\*AA, @prot_seqs);
}
if ($dna_only)
{
	&gjoseqlib::print_alignment_as_fasta(\@gene_seqs);
}
if ($aa_only)
{
	&gjoseqlib::print_alignment_as_fasta(\@prot_seqs);
}


chdir ($base);
unless ($keep_temp){system "rm -rf $tmp";}



##############sub call_non_pssm_features####################
#   Calles non pssm-based features that can be called based on 
#   protein coordinates.
#
#   call_non_pssm_features ($non_pssm_feat, $contigH)
#
#-----------------------------------------------------------
sub call_non_pssm_features
{
	my ($featH, $contigH) = @_;
	my @seq_data;
		
	foreach (keys %$featH)
	{
		my $feat          = $_;
		my $anno          = $featH->{$feat}->{ANNO};
		my $min           = $featH->{$feat}->{MIN};
		my $max           = $featH->{$feat}->{MAX};
		my $aa            = $featH->{$feat}->{AA};
		my $start_offset  = $featH->{$feat}->{START_OFFSET};
		my $stop_offset   = $featH->{$feat}->{STOP_OFFSET};
		my $stop_offset   = $featH->{$feat}->{STOP_OFFSET};
		my $feature_type  = $featH->{$feat}->{TYPE};
		my $symbol        = $featH->{$feat}->{SYMBOL};


		foreach (keys %{$featH->{$feat}->{COORD}})
		{
			my $contig = $_;
			
			if (($featH->{$feat}->{COORD}->{$contig}->{START}) && ($featH->{$feat}->{COORD}->{$contig}->{STOP}))
			{	
				my @starts = @{$featH->{$feat}->{COORD}->{$contig}->{START}};
				my @stops  = @{$featH->{$feat}->{COORD}->{$contig}->{STOP}};			
			
				for my $i (0..$#starts)
				{			
					my $start = $starts[$i]; 
					for my $j (0..$#stops)
					{
						my $stop = $stops[$j];
						if ($stop)
						{
							my ($strand, $begin, $end) = 0;
							if ($start < $stop)
							{
								$strand = "+";
								$begin = ($start += $start_offset);
								$end   = ($stop  -= $stop_offset);													
							}
							elsif ($start > $stop)
							{
								$strand = "-";
								$begin = ($start -= $start_offset);
								$end   = ($stop  += $stop_offset);						
							}
					
					
							if ((abs($begin - $end) <= $max) && (abs($begin - $end) >= $min))
							{
								print STDERR "\tCalling non-PSSM feature\t$anno\tbegin: $begin\tend: $end\n"; 
								my $nt = &gjoseqlib::DNA_subseq($contigH{$contig}, $begin, $end); 
								my $prot;
								if ($aa){ $prot = &gjoseqlib::translate_seq( $nt );}
								if ($prot =~ /\*(?!$)/)# It will not record a position-called protein with stops
								{
									print STDERR "\tInternal stop(s) found in non-pssm feature. No assignment made for $anno\n";
								}
								else
								{
									push @seq_data, ([$contig, $begin, $end, $anno, $strand, $feat, $nt, $prot, $feature_type, $symbol]);
								}
							}
						}
					}
				}
			}
		}
	}
	return @seq_data;
}













##########################sub crop_to_stop_codon###########
#
# Returns new gene boundaries, ignoring everything after the stop codon.
# Returns gene boundaries with stop codon by convention
#
#  $gene_end = crop_to_stop_codon($gene_start, $gene_end ,$hseq)
#  hseq is the AA sequence blast match containing the stop character (*)
#
#----------------------------------------------------------
sub crop_to_stop_codon
{
	my ($from, $to, $hseq) = @_; 		
	my $gene_end;
	$hseq =~ s/\*.+//g;
	my $len = length ($hseq); 
	my $to_keep = ($len * 3); 
	if ($from < $to){$gene_end = ($from + (($to_keep + 3) - 1))}
	if ($from > $to){$gene_end = (($from - ($to_keep) - 3) + 1) }
	return $gene_end;
}
###########################################################



##########################sub scan_to_stop_codon###########
# Looks beyond the end of the called blast boundaries to find the next stop codon
# returns the new gene coordinates, which includes the stop codon by convention
# The contig is the DNA string.  It is needed to ensure that you don't fall off of either end.
#
#   ($gene_begin, $gene_end)= scan_to_stop_codon ($to, $from, $contig);
#
#----------------------------------------------------------

sub scan_to_stop_codon
{
	my ($from, $to, $contig) = @_;
	my $len = length $contig; 
	my $end = $to;
		
	if ($from < $to)
	{		
		for (my $i = $end; $i <= ($len - 3); $i += 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i +1), ($i + 3) ); #<--- these are equivalent.
			my $aa = &gjoseqlib::translate_codon( $codon ); 
			
			if ($aa =~ /\*/) 
			{
				$end = ($i + 3);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			elsif ($aa =~ /x/i) 
			{
				print STDERR "\tC-term Extension X found was: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$end = ($i + 3);
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
			}		
		}
	}

	## Reverse:

	elsif ($from > $to)
	{		
		for (my $i = $end; $i > 3; $i -= 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i -1), $i - 3); #<--- these are equivalent.
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /\*/) 
			{
				$end = ($i - 3);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$aa\t$codon\n";

				last;	
			}		
			elsif ($aa =~ /x/i) 
			{
				print STDERR "\tC-term Extension X found was: $to\tnow: $end\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$end = ($i - 3);
				print STDERR "\tC-term Extended\twas: $to\tnow: $end\t$codon\t$aa\n";
			}		
		}
	}	
	return ($from, $end);
}
###########################################################


##########################sub scan_to_met_start#############
#
# Looks upstream for a new Met start codon.
# Returns the start position for the new Met codon.
# Does not consider alternative start codons.
#  ($new_gene_begin, $new_gene_end) = scan_to_met_start($gene_start, $gene_end, $contig)
# 
#-----------------------------------------------------------
sub scan_to_met_start
{
	my ($from, $to, $contig) = @_;
	my $len = length $contig; 
	my $start = $from;

	if ($from < $to)
	{		
		for (my $i = ($start - 3); $i >= 0; $i -= 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, $i, ($i + 2)); 
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /m/i) 
			{
				$start = $i;  #move the end position to the end of the stop codon and quit.
				print STDERR "\tN-term Extension Met found: was: $from\tnow: $start\t$aa\t$codon\n";
				last;	
			}		
			elsif ($aa =~ /x|\*/i) 
			{
				print STDERR "\tN-term Extension stopping was: $from\tnow: $start\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$start = $i;
				print STDERR "\tN-term Extended\twas: $from\tnow: $start\t$codon\t$aa\n";
			}		
		}
	}
	
	### Reverse:
	if ($from > $to)
	{		
		
		for (my $i = ($start + 1); $i <= ($len - 2); $i += 3) # This is zero indexed
		{
			my $codon = &gjoseqlib::DNA_subseq($contig, ($i + 2), $i); 
			my $aa = &gjoseqlib::translate_codon( $codon ); 

			if ($aa =~ /m/i) 
			{
				$start = ($i + 2);  #move the end position to the end of the stop codon and quit.
				print STDERR "\tN-term Extension Met found: was: $from\tnow: $start\t$aa\t$codon\n";
				last;	
			}		
			elsif ($aa =~ /x|\*/i) 
			{
				print STDERR "\tN-term Extension stopping was: $from\tnow: $start\t$codon\t$aa\n";
				last;	
			}		
			else
			{
				$start = ($i + 2);
				print STDERR "\tN-term Extended\twas: $from\tnow: $start\t$codon\t$aa\n";
			}		
		}
	}	
	return ($start, $to);
}
###########################################################



##########################sub matching_tblastn_hsps_json#####
# Find the best hsp from a tblastn blast result in json format
  
# Start with this:
#	use JSON::XS;
#   open (IN, "tblastn -outfmt 13 -db $s_file -in_pssm $pssm_file |");
#	my $blast = decode_json(scalar read_file(\*IN));	
#
#
#  (@results, $best_bit) = matching_tblastn_hsps_json($blast, $bitscore_cutoff, $coverage_cutoff);
#
#	where @results is a hash reference of every hsp meeting the bit score and coverage cutoff
#
#----------------------------------------------------------
sub matching_tblastn_hsps_json
{
	my ($blast, $bit_cutoff, $cov_cutoff) = @_;
	my $results = {};
	my @data;
	my $best_bit = 0;
	
	my @output = $blast->{BlastOutput2};
	for my $i (0..$#output)
	{
		my $iterations = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations};	
		my $nitr = scalar @$iterations;

		for my $j (0..($nitr -1 ))
		{
			if ($blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{iter_num} = 1)
			{
				my $hits = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits};
				my $nhits = scalar @$hits;
			
				for my $k (0..($nhits -1))
				{
					# okay for reference, if there are two contigs, the hsps will be split into two different "hits"
				
					# get the contig ID.  If i understand this correctly, the title has to be the same per HSP, so the hardcoded zero should be ok.
					my $contig = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{description}->[0]->{title};
					$contig =~ s/\s.+//g;								

					my $hsps = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps};
					my $nhsps = scalar @$hsps;
															
					for my $l (0..($nhsps -1))
					{
						$results->{contig}    = $contig;
						$results->{bit}       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{bit_score};						
						$results->{hseq}      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hseq};
						$results->{hit_from}  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_from};
						$results->{hit_to}    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_to};
						$results->{frame}     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{hit_frame};
						$results->{q_from}    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{query_from};
						$results->{q_to}      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{query_to};
						$results->{ali_len}   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{query_to};
						$results->{e_val}     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{evalue};
						$results->{hsp_num}   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps}->[$l]->{num};
				
						my $cov = ((length $results->{hseq})/(abs($results->{q_to} - $results->{q_from})+1));	
						$results->{cov} = $cov;

						if (($results->{bit} > $bit_cutoff) && ($cov > $cov_cutoff))
						{	
							push @data, $results;
						}
					
						if (($results->{bit} > $best_bit) && ($results->{bit} > $bit_cutoff))
						{
							$best_bit = $results->{bit};
						}
						$results = {};
					}
				}
			}
		}
	}
	return (\@data, $best_bit);
}
###########################################################



##########################sub get_blastn_bit_##############
# Dig the bit score out of the json format blastn
  
# Start with something like this:
#	use JSON::XS;
#		open (IN, "blastn -query $rep_file -subject $s_file -evalue 0.5 -reward 2 -penalty -3 -word_size 11 -outfmt 13 -soft_masking false |"); 
#		my $blastn = decode_json(scalar read_file(\*IN));	
#
#
#  $bit_score = get_blastn_bit($json);
#----------------------------------------------------------
sub get_blastn_bit 
{
    my $blast = shift @_;
    my $best_bit = 0;
    
    # Check if BlastOutput2 exists and is an array reference
    return 0 unless $blast->{BlastOutput2} && ref($blast->{BlastOutput2}) eq 'ARRAY';
    
    # Iterate through each BlastOutput2 element
    for my $output (@{$blast->{BlastOutput2}}) 
    {
        next unless $output->{report} && $output->{report}->{results};
        
        my $results = $output->{report}->{results};
        next unless $results->{bl2seq} && ref($results->{bl2seq}) eq 'ARRAY';
        
        # Iterate through bl2seq array
        for my $bl2seq_item (@{$results->{bl2seq}}) 
        {
            next unless $bl2seq_item->{hits} && ref($bl2seq_item->{hits}) eq 'ARRAY';
            next if @{$bl2seq_item->{hits}} == 0; # Skip empty hits arrays
            
            # Iterate through hits
            for my $hit (@{$bl2seq_item->{hits}}) 
            {
                next unless $hit->{hsps} && ref($hit->{hsps}) eq 'ARRAY';
                
                # Iterate through hsps
                for my $hsp (@{$hit->{hsps}}) 
                {
                    if ($hsp->{bit_score} && $hsp->{bit_score} > $best_bit) 
                    {
                        $best_bit = $hsp->{bit_score};
                    }
                }
            }
        }
    }
    
    return $best_bit;
}
###########################################################



