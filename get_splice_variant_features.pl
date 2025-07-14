#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use JSON::XS;
use File::Slurp;
use IPC::Run qw(run);
use Cwd;
use gjoseqlib;
use Getopt::Long::Descriptive;
use List::Util qw(max min);



my $program_description = <<'END_DESCRIPTION';
This program performs feature calling for spliced viral proteins.  It reads and writes GTO files.
It works by using a set of hand-curated nucleotide sequences.  The program first finds 
a BLASTn match to one of the reference sequqences, and then evaluates the match at the spice donor (SD)
and splice acceptor (SA) sites.  There must be conservation in the splice donor and splice acceptor sites
for it to work.

The reference sequences are formatted in their fasta headers as follows:

>valid sequence ID from a repo
*space* 
SD:SD_Region_Start-SD_Region_End;Last_nt_of_SD 
*space*
SA:SA_Region_Start-SA_Region_End;First_nt_of_SA 





END_DESCRIPTION

my $default_data_dir = $ENV{LOVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

my ($help, $tmp); 
my($opt, $usage) = describe_options(
    "\n$program_description\nUsage: %c %o",  
				    ["input|i=s"            => "Input GTO"],
				    ["output|o=s"           => "Output GTO"],
				    ["cov|c=i"              => "Overall Minimum BLASTn percent query coverage (D = 95)", { default => 95 }],
				    ["id|p=i"               => "Overall Minimum BLASTn percent identity  (D = 95)", { default => 95 }],
		            ["threads|a=i"          => "Threads for the BLASTN (D = 24))", { default => 24 }],
				    ["json|j=s"             => "Full path to the JSON opts file", {default => "$default_data_dir/Viral_PSSM.json"}],
				    ["dir|d=s"              => "Full path to the directory hand curated transcripts", {default => "$default_data_dir/Splice-Variants"}],
				    ["tmp|t=s"              => "Declare name for temp dir (D = randomly named in cwd)"], 
				    ["help|h"               => "Show this help message", { shortcircuit => 1 } ],
				    ["debug|b"              => "Enable debugging"],
				    ["seqs|s"             => "Dump sequences to STDERR"]
);


print $usage->text and exit if $opt->help;
die($usage->text) if @ARGV != 0;

if ($opt->tmp){$tmp = $opt->tmp;}
else {$tmp = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1))}
print STDERR "Tempdir=$tmp\n" if $opt->debug;

my $dir = $opt->dir;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";
$genome_in->{features}->[0] or die "No features in GTO\n"; 

my $base = getcwd;


# We read the GTO to get the family name
my %pssm_fam;
for my $i (0 .. $#{$genome_in->{features}}) 
{
	if (($genome_in->{features}->[$i]->{type} =~ /CDS/) && ($genome_in->{features}->[$i]->{family_assignments}->[0]->[3] =~ /LowVan/))
	{
		$pssm_fam{$genome_in->{features}->[$i]->{family_assignments}->[0]->[0]}++;
	}
}

die "More than one viral family of PSSMs in GTO\n" if scalar(keys %pssm_fam) > 1;
my $fam = (keys %pssm_fam)[0];
$fam or die "GTO has no annotations from LowVan Annotation tool\n"; 


# Next we read the JSON opts file to see if there are any spliced features that we need to find
my $json      = decode_json(read_file($opt->json));
$genome_in or die "Error reading json protein feature data";

my @to_analyze;
foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	if ($json->{$fam}->{features}->{$prot}->{special} eq "splice")
	{
		my $anno = $json->{$fam}->{features}->{$prot}->{anno};
		my $feature_type = $json->{$fam}->{features}->{$prot}->{feature_type};
		push @to_analyze, ([$prot, $anno, $feature_type]);
	}
}

# If there are spliced features then we set up the blast
print STDERR "\n\nSearching for spliced features\n------------------------------------\n"; 
print STDERR "\n$fam detected from GTO\n\n"; 

if (scalar @to_analyze)
{
	my %features;

	mkdir ($tmp); 
	chdir ($tmp);
	
	my $contigs = $genome_in->{id}."."."contigs"; 
	$genome_in->write_contigs_to_file($contigs);
	
	#make the blastn db in the temp dir.
	my $make_db = run("makeblastdb -dbtype nucl -in $contigs >/dev/null");

	if (!$make_db)
	{
  	 print STDERR "get_splice_variant_features:  makeblastdb failed with rc=$?. Stdout:\n";
	}
	
	# create the GTO analysis event.
	my $event = {
 	   tool_name => "LowVan Splice Variant Features",
  	   execution_time => scalar gettimeofday,
	};
	my $event_id = $genome_in->add_analysis_event($event);


	#cycle through the splice variant features and search for them one at a time with blastn.
	for my $i (0..$#to_analyze)
	{
		my $name = $to_analyze[$i][0];
		my $anno = $to_analyze[$i][1]; 
		my $ft   = $to_analyze[$i][2];
		
		print STDERR "Analyzing $name\n"; 
		my $query = "$dir/$fam/$name.fasta"; 
		run ("cp $query ."); 

		my $make_db2 = run("makeblastdb -dbtype nucl -in $name.fasta >/dev/null");

			if (!$make_db2)
			{
  				 print STDERR "get_transcript_edited_features:  makeblastdb failed with rc=$?. Stdout:\n";
			}

			my @blast_parms = (
		      "-query",         $query,
		      "-db",            $contigs, #### fix this in the original program
		      "-evalue",        0.5,
		      "-reward",          2,
		      "-penalty",        -3,
		      "-word_size",      28,
		      "-outfmt",          15,
		      "-soft_masking",   "false",
		      "-dust",           "no",
		      "-perc_identity",  $opt->id,   
              "-qcov_hsp_perc",  $opt->cov,
		      "-num_threads",    $opt->threads);
		      

		my $do_blast = run(["blastn", @blast_parms], ">", "$name.json", "2>", "$name.blastn.stderr.txt");
		open (IN, "<$name.json"), or warn "Cannot open JSON BLASTn output file $name.json\n";
		my $results = decode_json(scalar read_file(\*IN));	
		close IN;

		
		#Gather in the best match.
		my $matches = best_blastn_match_by_loc($results);  #removed id, cov, and gap thresholds from here  
		
		#Create a hash of valid Sequence Donor and Sequence Acceptor sites. 
		open (IN, "<$name.fasta"), or die "cannot open query sequence file to find valid SDs and SAs.\n";
		my @refs = &gjoseqlib::read_fasta(\*IN);
		my (%SDs, %SAs);
		for my $i (0..$#refs)
		{
			my ($sd_s, $sd_e, $sd_l, $sa_s, $sa_e, $sa_f) = $refs[$i][1] =~ /SD:(\d+)-(\d+);(\d+) SA:(\d+)-(\d+);(\d+)/;			
			my $SD = uc(substr($refs[$i][2], ($sd_s - 1), ($sd_e - $sd_s)+1));
			my $SA = uc(substr($refs[$i][2], ($sa_s - 1), ($sa_e - $sa_s)+1));
			$SDs{$SD} = 1;
			$SAs{$SA} = 1;
		}
		
		
		foreach (keys %$matches)
		{
			my $sid = $_;
			foreach (keys %{$matches->{$sid}})
			{
				my $from    = $_; 
				my $to      = $matches->{$sid}->{$from}->{TO};
				my $sseq    = $matches->{$sid}->{$from}->{HSEQ};
				my $qseq    = $matches->{$sid}->{$from}->{QSEQ};
				my $iden    = $matches->{$sid}->{$from}->{IDEN};
				my $ali_len = $matches->{$sid}->{$from}->{ALI_LEN};				
				my $qfrom   = $matches->{$sid}->{$from}->{QFROM};				
				my $qto     = $matches->{$sid}->{$from}->{QTO};				
				
				# Determine strand orientation
        		my $is_reverse = ($from > $to);
				my $strand = $is_reverse ? "-" : "+";

				my $pid     =  (($iden/$ali_len) * 100);				
				my $qcov    =  (($ali_len/(length $qseq)) * 100); 										
				my $qtitle  =  $matches->{$sid}->{$from}->{QTITLE};
				
				print STDERR "Best Match: $qtitle\n\n";
				
				#If all inclusion critreria are met:
				if ( ($pid >= $opt->id) && ($qcov >= $opt->cov) )
				{
					
					# Parse splice site coordinates from query title
					my ($sd_start, $sd_end, $sd_last, $sa_start, $sa_end, $sa_first) = $qtitle =~ /SD:(\d+)-(\d+);(\d+) SA:(\d+)-(\d+);(\d+)/;
					unless ($sd_start, $sd_end, $sd_last, $sa_start, $sa_end, $sa_first){warn "$name incorrectly formatted fasta header in query\n$qtitle\n"; }


					# Compute "relative" dna coordinates from the blastn alignment
					# These will always be in the (+) direction because the references must be in the plus.			
					my $rel_sd_start = $sd_start - $qfrom;
					my $rel_sd_end   = $sd_end - $qfrom;
					my $rel_sa_start = $sa_start - $qfrom;
					my $rel_sa_end   = $sa_end - $qfrom;
										
					# Ensure coordinates are within bounds and positive
					$rel_sd_start = max(0, $rel_sd_start);
					$rel_sd_end   = min($ali_len - 1, max(0, $rel_sd_end));
					$rel_sa_start = max(0, $rel_sa_start);
					$rel_sa_end   = min($ali_len - 1, max(0, $rel_sa_end));
					
					# Extract splice site aa sequences for validation
					#my $ref_SD = substr($qseq, $rel_sd_start, ($rel_sd_end - $rel_sd_start) + 1);
					#my $ref_SA = substr($qseq, $rel_sa_start, ($rel_sa_end - $rel_sa_start) + 1);
					
					my $sub_SD = substr($sseq, $rel_sd_start, ($rel_sd_end - $rel_sd_start) + 1);
					my $sub_SA = substr($sseq, $rel_sa_start, ($rel_sa_end - $rel_sa_start) + 1);
										
					
					if ((exists $SDs{$sub_SD}) && (exists $SAs{$sub_SA}))
					#if ((uc($ref_SD) eq uc($sub_SD)) && (uc($ref_SA) eq uc($sub_SA))) 
					{
						#print "Match found for SD and SA sites\n";
						if ($opt->seqs)
						{
							print STDERR "SD Match: $sub_SD\n";												
							print STDERR "SA Match: $sub_SA\n\n";
						}					
				
						# Calculate splice coordinates in the subject sequence
						my $splice_left_end    = $rel_sd_start + ($sd_last - $sd_start);
						my $splice_right_start = $rel_sa_start + ($sa_first - $sa_start);
						
						# Ensure coordinates are valid
						$splice_left_end    = max(0, min($ali_len - 1, $splice_left_end));
						$splice_right_start = max(0, min($ali_len - 1, $splice_right_start));
						
						# Create the splice variant by extracting left and right portions
						my $left_ungapped  = substr($sseq, 0, $splice_left_end + 1);
						my $right_ungapped = substr($sseq, $splice_right_start);
						
						# Remove gaps from the sequences
						$left_ungapped  =~ s/-//g;
						$right_ungapped =~ s/-//g;
						
						# I need to translate the left and right side to check them.
						my $trans_left = &gjoseqlib::translate_seq($left_ungapped);
						
			
						if ($trans_left =~ /\*/)
						{
							warn "Stop codon detected in 5-prime matching end of splice for: $name, skipping.\n";
							print STDERR "$left_ungapped\n$trans_left\n\n"; 
							#next; 
						}						
						
						# Create spliced sequence
						my $spliced_dna = $left_ungapped . $right_ungapped;
				
						# Translate sequences
						my $aa_splice = &gjoseqlib::translate_seq($spliced_dna);
				
						# Remove stop codons and downstream sequence
						$aa_splice =~ s/(\*)(.+)/$1/g;
						
						if ($opt->seqs)
						{	
							print STDERR "LEFT:\n$left_ungapped\n\nRIGHT:\n$right_ungapped\n\n"; 
							print STDERR "$aa_splice\n\n"; 
						}
						my $total_len = ((length ($aa_splice)) * 3);
						my $len_left  = length($left_ungapped); 
						my $len_right = length($right_ungapped); 
						my $cropped_len_right =  ($total_len - $len_left);					
						
						#Get genomic subject coordinates
						my (@left_tuple, @right_tuple, @loc);
						
						# Above I use "left" and "right" as terms that are relative to the
						# BLAST alignment. Since all references are in the "forward" 
						# orientation, this is fine (5'-3') for a subject on the forward strand.
						# In a reverse strand match, our notion of "left" is actually the 
						# right side and vice versa because the blast alignment is the RC. 
						# I correct for this when i build the tuple of location data for 
						# the two sections below.
						
						if ($is_reverse)
						{
							@left_tuple  = ([$sid, (($to + $len_right) - 1), $strand, $cropped_len_right]);
							@right_tuple = ([$sid, $from, $strand, $len_left]);						
							@loc = ([@right_tuple, @left_tuple]);
						}
						else
						{
							@left_tuple  = ([$sid, $from, $strand, $len_left]);
							@right_tuple = ([$sid, (($to - $len_right) + 1), $strand, $cropped_len_right]);
							@loc = ([@left_tuple, @right_tuple]);
						}
	
						#Debug locs.  Keeping this here for now.  
						#I have double checked this in the fwd, rev strands
						#Also double checked and is compatible with rast-export-genome.
						
						
						my $feature = {
							type        => $ft,
							contig      => $sid,
							aa_sequence => $aa_splice,
							location    => @loc,
							product     => $anno,
							pssm        => ([[$fam, $name, $anno, "LowVan Splice Variant Feature"]]),
						};
						push(@{$features{$ft}}, $feature)
	
					}
				}
			}
		}
	}

	if (%features)
	{
		foreach (keys %features)
		{
			my $type = $_; 
					
			foreach (@{$features{$type}})
			{
				my $data = $_;
	
				if ($type eq 'CDS' || $type eq 'mat_peptide')
				{
					my $p = {
						-id	                 => $genome_in->new_feature_id($type),
						-type 	             => $type,
						-location 	         => $data->{location},
						-analysis_event_id 	 => $event_id,
						-annotator           => 'LowVan Splice Variant Feature',
						-protein_translation => $data->{aa_sequence},
						-function            => $data->{product},
						-family_assignments  => $data->{pssm},
						};				
					$genome_in->add_feature($p);
				}
			}
		}			
		chdir ($base);
		$genome_in->destroy_to_file($opt->output);			
 	}
 	else
 	{
		#handle condition where there should have been a blast match, but none was found.
		print STDERR "Splice variant features should exist for: $fam, but none were found.\n";
		chdir($base);
		$genome_in->destroy_to_file($opt->output);			
	}
 }
 
 else 
 {
 	print STDERR "No splice variant features for: $fam\n";
 	chdir($base);
 	$genome_in->destroy_to_file($opt->output);			
 }






















##########################sub best_blastn_match_by_loc##############
# Find the best blast match(s) for a blastn json
#
#   This will return the best blastn matched sequence by location.  Only one HSP within 
#   A distance defined as > abs(alignment length) is returned per contig.  
#   (i.e., you can have more than one "best" blast match per contig, but they
#   can't overlap.  The amount of overlap could be turned into a parameter, 
#   but i just used the ali length, which seemed reasonable. 
#
#   NOTE: this reads the -db formatted json output not the -subject [fasta] formatted version. 
#   The returned JSONs are slightly different 
#
#   Relevant paramaters such as %id %coverage should be set in the blast command line options,
#   or post-processed.
#
#   usage:
#   $hash = best_blastn_match_by_loc($blastn_json);
# 
#   The returned hash reference is in the following format:
#     
#    SubjectID->{HitFrom}->{TO}->{HitTo} 
#             ->{HitFrom}->{HSEQ}->{SubjectSeq}
#             ->{HitFrom}->{QSEQ}->{QuerySeq}
#             ->{HitFrom}->{BIT}->{BitScore} 
#             ->{HitFrom}->{IDEN}->{NumIdentities}         ### need to add this
#             ->{HitFrom}->{ALI_LEN}->{Alignment_Length}   ### need to add this
#             ->{HitFrom}->{GAPS}->{gaps}
#
# 
#----------------------------------------------------------
sub best_blastn_match_by_loc
{
	my ($blast)  =  @_;
	my $matches = {};
	
	my @output = @{$blast->{BlastOutput2}};
	for my $i (0..$#output)
	{
		my @hits = @{$blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}};
		for my $j (0..$#hits)
		{
			my @hsps = @{$blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}};
			for my $k (0..$#hsps)
			{
				my $ident     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{identity};
				my $ali_len   = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{align_len};
				my $qlen      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{query_len};
				my $qtitle    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{query_title};
				my $gaps      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{gaps};
				my $sid       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{description}->[$k]->{title};
				my $hit_from  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_from};
				my $bit       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{bit_score};
				my $hit_to    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_to};
				my $qseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{qseq};
				my $hseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hseq};
				my $qfrom     = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{query_from};
				my $qto       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{query_to};

								
				my ($pid, $qcov);
				if ($ident && $ali_len && $qlen) # this ensures that we got search results.
				{
					# if we have already seen a match in this contig:
					if (exists $matches->{$sid}) 
					{							
						foreach (keys %{$matches->{$sid}})
						{
							my $loc = $_; #hit_from location that was seen previously
							
							#if the match is in roughly in the same location as seen before, and it has a better bitscore, 
							#then they are considered the same match.  We delete the orignal and keep the one with the better bit score.							
							if ((abs($hit_from - $loc) < $ali_len) && ($bit > $matches->{$sid}->{$loc}->{BIT}))
							{
								delete $matches->{$sid}->{$loc};
								$matches->{$sid}->{$hit_from}->{BIT}      = $bit;
								$matches->{$sid}->{$hit_from}->{TO}       = $hit_to;
								$matches->{$sid}->{$hit_from}->{QSEQ}     = $qseq;
								$matches->{$sid}->{$hit_from}->{HSEQ}     = $hseq;
								$matches->{$sid}->{$hit_from}->{IDEN}     = $ident;
								$matches->{$sid}->{$hit_from}->{ALI_LEN}  = $ali_len;
								$matches->{$sid}->{$hit_from}->{GAPS}     = $gaps;
								$matches->{$sid}->{$hit_from}->{QTITLE}   = $qtitle;
								$matches->{$sid}->{$hit_from}->{QTO}      = $qto;
								$matches->{$sid}->{$hit_from}->{QFROM}    = $qfrom;
							
							}
						}
					}
					#if we have not seen a match in this location before add it:
					else  
					{
						$matches->{$sid}->{$hit_from}->{BIT}      = $bit;
						$matches->{$sid}->{$hit_from}->{TO}       = $hit_to;
						$matches->{$sid}->{$hit_from}->{QSEQ}     = $qseq;
						$matches->{$sid}->{$hit_from}->{HSEQ}     = $hseq;
						$matches->{$sid}->{$hit_from}->{IDEN}     = $ident;
						$matches->{$sid}->{$hit_from}->{ALI_LEN}  = $ali_len;
						$matches->{$sid}->{$hit_from}->{GAPS}     = $gaps;
						$matches->{$sid}->{$hit_from}->{QTO}      = $qto;
						$matches->{$sid}->{$hit_from}->{QFROM}    = $qfrom;
						$matches->{$sid}->{$hit_from}->{QTITLE}   = $qtitle;

					}				
				}
			}
		}
	}
	return $matches;
}

###########################################################













