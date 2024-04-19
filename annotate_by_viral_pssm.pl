#! /usr/bin/env perl

use strict;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;
use Cwd;
use gjoseqlib;

my $usage = 'annotate_by_viral_pssm.pl [options] -i subject_contig(s).fasta 

		-h help
		-i input subject contigs in fasta format
		-t declare a temp file (d = random)
		-tax declare a taxonomy id (D = 11158 )
		-g Genome name (D = Paramyxoviridae);
		
		-min minimum contig	length (d = 1000) # otherwise the genome is rejected
		-max maximum contig length (d = 30000) # for reference Measles is 15894 and beilong is 19,212

        -opt Options file in JSON format which carries data for match (D = /home/jjdavis/Viral_PSSM.json)
		-l Representative contigs directory (D = /home/jjdavis/bin/Viral-Rep-Contigs)
		-p Base directory of PSSMs   (D = /home/jjdavis/bin/Viral-PSSMs)
	      Note that this is set up as a directory of pssms
	      right now this is hardcoded as: "virus".pssms within this directory.

	    -ks if a pssm match extends beyond a stop codon it will keep the entire region of the match including the stop. 
	      
	   Debug parms (turns off output file generation)
	   	-tmp keep temp dir
		-no no output files generated, to be used in conjunction with one of the following:
		    -dna print only genes to STDOUT 
		    -aa print proteins to STDOUT
		    -tbl print only feature table to STDOUT
		    -ctbl [file name] concatenate table results to a file (for use with many genomes)

';

my ($help, $opt_file, $contig_file, $tmp, $tax, $keep_stop, $genome_name, $cdir, $pdir, $keep_temp, $min_len, $max_len, $aa_only, $dna_only, $tbl_only, $no_out, $ctbl);

my $opts = GetOptions( 'h'         => \$help,
                       'tmp'       => \$keep_temp,
                       'no'        => \$no_out,
                       'aa'        => \$aa_only,
                       'dna'       => \$dna_only,
                       'tbl'       => \$tbl_only,
                       'ctbl=s'    => \$ctbl,
                       'i=s'       => \$contig_file,
                       'tax=s'     => \$tax,
                       't=s'       => \$tmp,
                       'g=s'       => \$genome_name,
                       'l=s'       => \$cdir,
                       'p=s'       => \$pdir,
                       'min=s'     => \$min_len,
                       'max=s'     => \$max_len,
                       'opt=s'     => \$opt_file,
                       'ks'        => \$keep_stop); 

if ($help){die "$usage\n";}
unless ($contig_file ){die "must declare an input subject file with -i \n\n$usage\n";}

# Name the Temp file:
# generates a random 20 digit string of 0-9a-f
unless ($tmp){$tmp .= sprintf("%x", rand 16) for 1..20;}
unless ($min_len){$min_len = 1000; }
unless ($max_len){$max_len = 30000; }
unless ($tax){$tax = "11158"; }
#unless ($name){$name = "Paramyxoviridae"; }
unless ($cdir){$cdir = "/home/jjdavis/bin/Viral-Rep-Contigs"; }
unless ($pdir){$pdir = "/home/jjdavis/bin/Viral-PSSMs"; }
unless ($opt_file){$opt_file = "/home/jjdavis/bin/Viral_PSSM.json"; }


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
	if ($genome_name)
	{
		open (AA,  ">$genome_name.faa");
		open (DNA, ">$genome_name.ffn");
		open (TBL, ">$genome_name.feature.tbl");
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
foreach (@reps)
{
	my $rep = $_;
	my $virus = $_;
	$virus =~ s/\..+//g; 
	
	my $rep_file = "$cdir/$rep";
	open (IN, "blastn -query $rep_file -subject $s_file -evalue 0.5 -reward 2 -penalty -3 -word_size 11 -outfmt 15 -soft_masking false |"); 
	my $blastn = decode_json(scalar read_file(\*IN));	
	close IN;

	my $match_bit = get_blastn_bit($blastn);  
	unless ($match_bit < $best_contig_bit)
	{
		$best_contig_bit = $match_bit;
		$best_virus_match = $virus;
	}
}	


my $virus = $best_virus_match;	
print STDERR "-----------------------\nMatching $virus\tBit = $best_contig_bit\n-----------------------\n"; 
opendir (DIR, "$pdir/$virus.pssms");
my @pssm_dirs = grep{$_ !~ /^\./}readdir(DIR);  # reads the directory of PSSM dirs
close DIR;

my @all_seqs;
my %positions;
#my ($upstream_ext, $downstream_ext, $keep_stop, $bit_cutoff, $cov_cutoff);

my $join = {};
my $non_pssm_feat = {};
foreach (@pssm_dirs)  #Each PSSM dir contains one or more PSSMs for a given homolog
{
	my $pssmdir = $_; 
	opendir (DIR, "$pdir/$virus.pssms/$pssmdir");
	my @pssm_files = grep{$_ !~ /^\./}readdir(DIR); # gets the pssms
	close DIR;

	# retrieve the protein specific options for the corresponding PSSM.
	unless ($options->{$virus}->{$pssmdir})
	{
		print STDERR "No data in JSON file for VIRUS: $virus\tPSSM: $pssmdir\n"; 
		next;
	}
	my $upstream_ext   = $options->{$virus}->{$pssmdir}->{upstream_ext};
	my $downstream_ext = $options->{$virus}->{$pssmdir}->{downstream_ext};
	my $bit_cutoff     = $options->{$virus}->{$pssmdir}->{bit_cutoff};
	my $cov_cutoff     = $options->{$virus}->{$pssmdir}->{coverage_cutoff};
	my $start_to_met   = $options->{$virus}->{$pssmdir}->{start_to_met};
		
	print STDERR "\t$virus\t$pssmdir\tbit\t$bit_cutoff\tcov\t$cov_cutoff\tkeep_stop\t$keep_stop\tupstream_ext\t$upstream_ext\tdownstream_ext\t$downstream_ext\n"; 		

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
		open (IN, "tblastn -outfmt 15 -db $s_file -in_pssm $pssm_file |");
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
	print STDERR "\tChoosing $best_pssm\t$best_bit\n"; 

	my $nhsps = scalar @$best_results;
	for my $i (0..($nhsps -1))
	{		
		my $contig = $best_results->[$i]->{contig};
		$contig =~ s/\s.+//g;								
		my $hseq = $best_results->[$i]->{hseq};
		$hseq =~ s/\-//g; # eliminate alignment gaps in hit seq
			
		print STDERR "\t$contig\t$best_results->[$i]->{anno}\t$best_results->[$i]->{hit_from}\t$best_results->[$i]->{hit_to}\t$best_results->[$i]->{frame}\t$best_results->[$i]->{bit}\n";
	
		
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
		#print STDERR "Original Coords: $from\tto\t$to\n"; 	

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
			($gene_begin, $gene_end)= scan_to_stop_codon ($to, $from, $contigH{$contig});
			print STDERR "\tScanning to stop codon Original Coords: $from\tto\t$to\tNew Coords: $gene_begin\tto\t$gene_end\n"; 	

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
			print STDERR "Scanning for Met start\n";
			($new_gene_begin, $new_gene_end) = scan_to_met_start($gene_begin, $gene_end, $contigH{$contig});
		}
		if ($new_gene_begin){$gene_begin = $new_gene_begin};	


		my $gene = &gjoseqlib::DNA_subseq($contigH{$contig}, $gene_begin, $gene_end );	
		my $protein = &gjoseqlib::translate_seq( $gene );

		if ($start_to_met){$protein =~ s/^[A-Z]/m/i;}
		
		# Set up the Paramyxo Join (if necessary)
		if (exists $options->{$virus}->{$pssmdir}->{paramyxo_join})
		{			
			$join->{$contig}->{$pssmdir}->{START}    = $gene_begin;
			$join->{$contig}->{$pssmdir}->{STOP}     = $gene_end;
			$join->{$contig}->{$pssmdir}->{PARTNER}  = $options->{$virus}->{$pssmdir}->{join_partner};
			$join->{$contig}->{$pssmdir}->{ORDER}    = $options->{$virus}->{$pssmdir}->{paramyxo_join};
			$join->{$contig}->{$pssmdir}->{ANNO}     = $options->{$virus}->{$pssmdir}->{new_anno};	
			$join->{$contig}->{$pssmdir}->{INSERT}   = $options->{$virus}->{$pssmdir}->{paramyxo_insert};	
			$join->{$contig}->{$pssmdir}->{PSSM}     = $best_pssm;
		}
		
	
	
		# Set up calling non-pssm features that are anchored to pssm coordinates
		if (exists $options->{$virus}->{$pssmdir}->{non_pssm_partner})
		{
			my @non_pssms = @{$options->{$virus}->{$pssmdir}->{non_pssm_partner}};
			for my $i (0..$#non_pssms)
			{
				my $start_coord;
				my $stop_coord;
				my $contig = $best_results->[$i]->{contig};
				my $feat = $non_pssms[$i];
				#if the current pssm match feature corresponds with the non-pssm start site
				if ($options->{$virus}->{$feat}->{begin}->{begin_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{START_OFFSET} = $options->{$virus}->{$feat}->{begin}->{begin_offset};
					if ($options->{$virus}->{$feat}->{begin}->{begin_pssm_loc} =~ /START/)
					{	
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_begin;
					}
					elsif ($options->{$virus}->{$feat}->{begin}->{begin_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{START}}, $gene_end;
					}
				}
				elsif ($options->{$virus}->{$feat}->{end}->{end_pssm} eq $pssmdir) 
				{
					$non_pssm_feat->{$feat}->{STOP_OFFSET} = $options->{$virus}->{$feat}->{end}->{end_offset};
					if ($options->{$virus}->{$feat}->{end}->{end_pssm_loc} =~ /START/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_begin;
					}
					elsif ($options->{$virus}->{$feat}->{end}->{end_pssm_loc} =~ /STOP/)
					{
						push @{$non_pssm_feat->{$feat}->{COORD}->{$contig}->{STOP}}, $gene_end;
					}
				}			
				$non_pssm_feat->{$feat}->{ANNO} = $options->{$virus}->{$feat}->{anno};
				$non_pssm_feat->{$feat}->{MIN}  = $options->{$virus}->{$feat}->{min_len};
				$non_pssm_feat->{$feat}->{MAX}  = $options->{$virus}->{$feat}->{max_len};			
				$non_pssm_feat->{$feat}->{AA}   = $options->{$virus}->{$feat}->{translate};
			}
		}
		
		#return matching sequence as a tuple.
		unless ($options->{$virus}->{$pssmdir}->{paramyxo_join} == 2)# the ORF2 sequence isn't complete and has to be merged after all of the other proteins have been found.
		{
			push @all_seqs, ([$best_results->[$i]->{contig}, $gene_begin, $gene_end, $best_results->[$i]->{anno}, $strand, $best_pssm, $gene, $protein]); 
		}	
	}
	print STDERR "-----------------------\n"; 
}


#create any joins if they exist:
if ($join)
{
	#this is kind of ugly because its taking a hash_ref for join and a hash for the contigs 
	my @tuples = join_orfs ($join, %contigH);
	push @all_seqs, @tuples;
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
		my $prot_id = "jim\|$tax\.$version\.CDS\.$count"; 
		
		#Trimming the star from the end of the protein sequence. 
		#This could be removed if its unwanted.
		my $protein = $_->[7];
		$protein =~ s/\*$//g;

		push @gene_seqs, ([$prot_id, $_->[3], $_->[6]]);
		push @prot_seqs, ([$prot_id, $_->[3], $protein]);
		
		# Original Feature Table looks like this:
		# genome_id	genome_name	accession	annotation	feature_type	patric_id	refseq_locus_tag	start	end	strand	na_length	gene	product	plfam_id	pgfam_id
		# For now, I will keep:
		# genome_id	genome_name	accession(contig)	annotation	feature_type	patric_id	start	end	strand	na_length	gene	product	plfam_id	pgfam_id

		my $na_len = length $_->[6];
		unless($no_out)
		{
			print TBL "$tax\.$version\t$genome_name\t$_->[0]\tJIM\tCDS\t$prot_id\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$_->[5]\t$_->[3]\n";
		}
		if($tbl_only)
		{
			print "$tax\.$version\t$genome_name\t$_->[0]\tJIM\tCDS\t$prot_id\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$_->[5]\t$_->[3]\n";
		}
		if($ctbl)
		{
			print CTBL "$tax\.$version\t$genome_name\t$_->[0]\tJIM\tCDS\t$prot_id\t$_->[1]\t$_->[2]\t$_->[4]\t$na_len\t$_->[5]\t$_->[3]\n";
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
		my $feat         = $_;
		my $anno         = $featH->{$feat}->{ANNO};
		my $min          = $featH->{$feat}->{MIN};
		my $max          = $featH->{$feat}->{MAX};
		my $aa           = $featH->{$feat}->{AA};
		my $start_offset = $featH->{$feat}->{START_OFFSET};
		my $stop_offset  = $featH->{$feat}->{STOP_OFFSET};

		print "####STARTOFFSET = $start_offset\n###STOP_OFFSET = $stop_offset\n\n";
		
		foreach (keys %{$featH->{$feat}->{COORD}})
		{
			my $contig = $_;
			
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
						my ($strand, $begin, $end);
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
					
						print "######STRAND = $strand\n Begin = $begin\n END = $end\n\n"; 
					
						if ((abs($begin - $end) <= $max) && (abs($begin - $end) >= $min))
						{
							my $nt = &gjoseqlib::DNA_subseq($contigH{$contig}, $begin, $end); 
							my $prot;
							if ($aa){ $prot = &gjoseqlib::translate_seq( $nt );}
							push @seq_data, ([$contig, $begin, $end, $anno, $strand, $feat, $nt, $prot]);
						}
					}
				}
			}
		}
	}
	return @seq_data;
}
############################################################










##########################sub join_orfs######################
# Takes the join hashref from the main loop, and the contigs Hash and joins the orfs
#    (usually orfs within the phosphoprotien)
# 
#     my @seq_data = join_orfs($join, %contigsH);
# 
# Returns an array of arrays that can be pushed into @all_seqs:
# [contig, start, end, anno, strand, matching_pssms, gene, protein]
# 
# The join will look like this:  [PotA start] |------>    [ProtA end]
#              		             [ORF-2 start]    |-----> [ORF-2 End]			
# 
# Where we need the first part of Protein A, spliced with ORF-2, to make full-length protein B.
# The PSSM for protein B is not a full-length protein because there is a frame jump.
# The part after the frame jump is what I am calling ORF-2 and it has a separate PSSM. 
#
#  NOTE that the returned "gene" is the unaltered DNA strand 
#  from the start of protA to the end of ORF-2.  Because Protein B is the result of 
#  "RNA editing" the inserted Gs from the stuttering RdRp will not appear in the sequence.
#
#  The protein, however, if spliced correctly will be the correct AA seq.
#-----------------------------------------------------------
sub join_orfs
{
	my ($join, $contigH) = @_;
	my @seq_data;
		
	foreach (keys %$join)
	{
		my $contig = $_;
		foreach (keys %{$join->{$contig}})
		{
			my $name = $_;
			if ($join->{$contig}->{$name}->{ORDER} == 1)	
			{
				my @partners = @{$join->{$contig}->{$name}->{PARTNER}};
				for my $i (0..$#partners)
				{
					my $orf2       = $join->{$contig}->{$name}->{PARTNER}->[$i];
					my $orf2_start = $join->{$contig}->{$orf2}->{START}; 
					my $orf2_end   = $join->{$contig}->{$orf2}->{STOP}; 
					my $orf1_start = $join->{$contig}->{$name}->{START}; 
					my $orf1_end   = $join->{$contig}->{$name}->{STOP}; 
					my $insert     = $join->{$contig}->{$orf2}->{INSERT}; 
					my $pssms      = "$join->{$contig}->{$name}->{PSSM}"."_JOIN_"."$join->{$contig}->{$orf2}->{PSSM}";
					
					#print Dumper $join;
					
					
					if (($orf1_start) && ($orf2_start) && ($orf1_end) && ($orf2_end))
					{	
						#Forward Strand:
						if ($orf1_end > $orf1_start)
						{
							#make sure the second orf start is within the first orf. 
							if (( $orf2_start > $orf1_start) && ($orf2_start <= $orf1_end))
							{
								my $gene = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf1_start, $orf2_end);			
								# make the spliced protein.
								
								my $ntA = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf1_start, $orf2_start); 
								my $ntB = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf2_start, $orf2_end);
							
								my $transcript = $ntA.$ntB; 
							
								#print STDERR "$ntA\n$insert\n$ntB\n"; 
								if ($insert)
								{
									print STDERR "ADDING INSERT\t$insert\n";
									$transcript = $ntA.$insert.$ntB;
								}
							
								my $protein = &gjoseqlib::translate_seq( $transcript );
								#print STDERR "####PROT: $protein\n"; 
								push @seq_data, ([$contig, $orf1_start, $orf2_end, $join->{$contig}->{$orf2}->{ANNO}, "+", $pssms, $gene, $protein]);
							}
							else
							{
								print STDERR "$name\t$orf2\tstart and stop do not overlap in the join\n"; 
							}
						}		
						#Reverse Strand:
						elsif ($orf1_end < $orf1_start)
						{
							#make sure the second orf start is within the first orf. 
							if (( $orf2_start < $orf1_start) && ($orf2_start >= $orf1_end))
							{
								my $gene = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf1_start, $orf2_end);					
								# make the spliced protein.
								my $ntA = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf1_start, $orf2_start);##<-- this splice
								my $ntB = &gjoseqlib::DNA_subseq($contigH{$contig}, $orf2_start, $orf2_end);
														
								my $transcript = $ntA.$ntB; 
							
								#print STDERR "$ntA\n$insert\n$ntB\n"; 
								if ($insert)
								{
									print STDERR "ADDING INSERT\t$insert\n";
									$transcript = $ntA.$insert.$ntB;
								}
							
								my $protein = &gjoseqlib::translate_seq( $transcript );

								my $protA = &gjoseqlib::translate_seq( $ntA );
								my $protB = &gjoseqlib::translate_seq( $ntB );
								my $protein = "$protA.$protB";
				
								push @seq_data, ([$contig, $orf1_start, $orf2_end, $join->{$contig}->{$orf2}->{ANNO}, "-", $pssms, $gene, $protein]);

							}
							else
							{
								print STDERR "$name\t$orf2\tstart and stop do not overlap in the join\n"; 
							}
						}
					}
				}
			}
		}
	}
	return @seq_data;
}
############################################################







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
	my ($begin, $end, $contig) = @_;
	my ($end_c, $begin_c);

	if ($begin < $end){$end_c = ($begin - 1)}
	elsif ($begin > $end){$end_c = ($begin + 1)}
	
	my $len = length($contig);
	while (($end_c > 0) && ($end_c <= $len)) 
	{
		if ($begin < $end){ $begin_c = ($end_c - 2);} 
		elsif ($begin > $end){ $begin_c = ($end_c + 2);} 
	
		my $codon = &gjoseqlib::DNA_subseq($contig, $begin_c, $end_c );	
		my $aa = lc(&gjoseqlib::translate_codon( $codon ));
		
		last if ($aa =~ /M/i);
		
		if (($aa =~ /\*/)|| ($aa =~ /x/i))#have to revert to the prior codon if you hit a stop or ambiguous codon.
		{
			if ($begin < $end){$begin_c += 3};
			if ($begin > $end){$begin_c -= 3};
			last;
		}
		if ($begin < $end){$end_c = ($begin_c - 1 )};
		if ($begin > $end){$end_c = ($begin_c + 1 )};
		print STDERR "\tN-term Extended:\tori_start:$begin\tnew_start:$begin_c\t$codon\t$aa\n";

	}							
	return ($begin_c, $end);
}
###########################################################



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
	my ($to, $from, $contig) = @_;

	my $end;
	if ($from < $to) {$end = ($to + 1);}
	elsif ($from > $to) {$end = ($to - 1);}

	my $contiglen = length($contig);
	
	#print STDERR "END=$end\t$contiglen\n";

	while (($end <= $contiglen) && ($end > 0))
	{
		my $next_c;
		if ($from < $to){$next_c = ($end + 2)}
		if ($from > $to){$next_c = ($end - 2)}
		
		my $codon = &gjoseqlib::DNA_subseq($contig, $end, $next_c );	
		my $aa = lc(&gjoseqlib::translate_codon( $codon ));
		
		if ($aa =~ /x/i) #revert to the previous codon if the current codon is ambiguous.
		{
			if ($from < $to){ $end -= 3;}
			if ($from > $to){ $end += 3;}
			last;
		}
		
		$end = $next_c;
		print STDERR "\tC-term Extended:was:$to\tnow:$next_c\t$codon\t$aa\n";
		
		last if ($aa =~ /\*/);    #if its not a stop codon or ambiguous, increment the end position.
		
		if ($from < $to){ $end = ($next_c + 1)}
		elsif ($from > $to){ $end = ($next_c - 1)}
	}
	return ($from, $end);
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
					#right now, I'm using the PSSM title as the annotation.
					my $anno = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{query_title};  

					my $hsps = $blast->{BlastOutput2}->[$i]->{report}->{results}->{iterations}->[$j]->{search}->{hits}->[$k]->{hsps};
					my $nhsps = scalar @$hsps;
															
					for my $l (0..($nhsps -1))
					{
						$results->{contig}    = $contig;
						$results->{anno}      = $anno;
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
	my $bit;
	my $best_bit = 0;

	#my @results = $blast->{BlastOutput2}->[0]->{report}->{results}->{bl2seq};	
	my @output = $blast->{BlastOutput2};
	for my $i (0..$#output)
	{
		my @results = $blast->{BlastOutput2}->[$i]->{report}->{results}->{bl2seq};

		for my $j (0..$#results)
		{
			my @hits = $blast->{BlastOutput2}->[$i]->{report}->{results}->{bl2seq}->[$j]->{hits};
			for my $k (0..$#hits)
			{
				my @hsps = $blast->{BlastOutput2}->[$i]->{report}->{results}->{bl2seq}->[$j]->{hits}->[$k]->{hsps}; 
				for my $l (0..$#hsps)
				{
					$bit = $blast->{BlastOutput2}->[$i]->{report}->{results}->{bl2seq}->[$j]->{hits}->[$k]->{hsps}->[$l]->{bit_score};
					if ($bit > $best_bit)
					{
						$best_bit = $bit;
					}
				}
			}
		}
	}
	return $best_bit;
}
###########################################################









