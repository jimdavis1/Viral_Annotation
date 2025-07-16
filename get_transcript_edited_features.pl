#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
#use P3DataAPI;
use JSON::XS;
use File::Slurp;
#use File::Path 'remove_tree';
use IPC::Run qw(run);
use Cwd;
use gjoseqlib;
use Getopt::Long::Descriptive;


my $program_description = <<'END_DESCRIPTION';
This program performs feature calling for transcript edited proteins.  It reads and writes GTO files.
It works by using a set of hand-curated transcripts in --dir as the queries.  --id, --gaps, and --cov refer 
to the strict inclusion criteria for enabling the mapping the nucleotides from the closest hand-curated transcript
onto the subject sequence.  When we enoucnter a blast match that is good, [defined by --lower_pid, --lower_pcov, and --eval], 
but not good enough to carry over the transcript-edited seqeunce, we call a partial_cds feature and the annotation becomes:
[Uncorrected . annotation string . encoding region]. It is considered a partial_cds feature because the translation would be 
be interrupted where the frame jump occurs, or shortly thereafter.

END_DESCRIPTION


my $default_data_dir = $ENV{LOVAN_DATA_DIR} // "/home/jjdavis/bin/Viral_Annotation";

my ($help, $tmp); 
my($opt, $usage) = describe_options(
    "\n$program_description\nUsage: %c %o",  
				    ["input|i=s"            => "Input GTO"],
				    ["output|o=s"           => "Output GTO"],
				    ["cov|c=i"              => "Minimum BLASTn percent query coverage (D = 95)", { default => 95 }],
				    ["id|p=i"               => "Minimum BLASTn percent identity  (D = 95)", { default => 95 }],
				    ["gaps|g=i"             => "Maximum number of allowable gaps (D = 2)", { default => 2 }],
				    ["e_val|e=f"            => "Maximum BLASTn evalue for considering any HSP (D = 0.5)", { default => 0.5 }],
				    ["lower_pid|lpi=i"      => "Lower percent identity threshold for a feature call without transcript editing correction (D = 80)", {default => 80}],
				    ["lower_pcov|lpc=i"     => "Lower percent query coverage for for a feature call without transcript editing correction (D = 80)", {default => 80}],
				    ["threads|a=i"          => "Threads for the BLASTN (D = 24))", { default => 24 }],
				    ["json|j=s"             => "Full path to the JSON opts file", {default => "$default_data_dir/Viral_PSSM.json"}],
				    ["dir|d=s"              => "Full path to the directory hand curated transcripts", {default => "$default_data_dir/Transcript-Editing"}],
				    ["tmp|t=s"              => "Declare name for temp dir (D = randomly named in cwd)"], 
				    ["help|h"               => "Show this help message", { shortcircuit => 1 } ],
				    ["debug|b"              => "Enable debugging"],
);


print $usage->text and exit if $opt->help;
die($usage->text) if @ARGV != 0;

if ($opt->tmp){$tmp = $opt->tmp;}
#else{$tmp .= sprintf("%x", rand 16) for 1..20;}
else {$tmp = File::Temp->newdir(CLEANUP => ($opt->debug ? 0 : 1))}
print STDERR "Tempdir=$tmp\n" if $opt->debug;

my $dir = $opt->dir;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";
$genome_in->{features}->[0] or die "No features in GTO\n"; 

my $base = getcwd;


# We read the GTO to get the family 
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


# Next we read the JSON to see if there are any transcript edited features that we need to find
my $json      = decode_json(read_file($opt->json));
$genome_in or die "Error reading json protein feature data";

my @to_analyze;
foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	if ($json->{$fam}->{features}->{$prot}->{special} eq "transcript_edit")
	{
		my $anno = $json->{$fam}->{features}->{$prot}->{anno};
		my $symbol = $json->{$fam}->{features}->{$prot}->{gene_symbol};
		my $feature_type = $json->{$fam}->{features}->{$prot}->{feature_type};
		push @to_analyze, ([$prot, $anno, $feature_type, $symbol]);
	}
}

# If there are transcript edited features then we set up the blast
if (scalar @to_analyze)
{
	my %features;
	print STDERR "\n\nSearching for transcript-edited features\n------------------------------------\n"; 
	print STDERR "\n$fam detected from GTO\n\n"; 

	mkdir ($tmp); 
	chdir ($tmp);
	
	my $contigs = $genome_in->{id}."."."contigs"; 
	$genome_in->write_contigs_to_file($contigs);
	
	#make the blastn db in the temp dir.
	my $make_db = run("makeblastdb -dbtype nucl -in $contigs >/dev/null");

	if (!$make_db)
	{
  	 print STDERR "get_transcript_edited_features:  makeblastdb failed with rc=$?. Stdout:\n";
	}
	
	# create the GTO analysis event.
	my $event = {
 	   tool_name => "LowVan Transcript Edited Features",
  	   execution_time => scalar gettimeofday,
	};
	my $event_id = $genome_in->add_analysis_event($event);

	
	
	#cycle through the transcript edited features and search for them one at a time with blastn.
	for my $i (0..$#to_analyze)
	{
		my $name   = $to_analyze[$i][0];
		my $anno   = $to_analyze[$i][1]; 
		my $ft     = $to_analyze[$i][2];
		my $symbol = $to_analyze[$i][3];

		
		print STDERR "\tAnalyzing $name\n\n"; 
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
		      "-evalue",        $opt->e_val,
		      "-reward",          2,
		      "-penalty",        -3,
		      "-word_size",      28,
		      "-outfmt",         15,
		      "-soft_masking",   "false",
		      "-dust",           "no",
		      "-perc_identity",  $opt->lower_pid,   
              "-qcov_hsp_perc",  $opt->lower_pcov,
		      "-num_threads",    $opt->threads);

		my $do_blast = run(["blastn", @blast_parms], ">", "$name.json", "2>", "$name.blastn.stderr.txt");
		
		open (IN, "<$name.json"), or warn "Cannot open JSON BLASTn output file $name.json\n";
		my $results = decode_json(scalar read_file(\*IN));	
		close IN;
			
				
		#Gather in the best match.
		my $matches = best_blastn_match_by_loc($results);  #removed id, cov, and gap thresholds from here  
		

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
				my $dashes  = ($sseq =~ tr/-//); 
				my $runs    = (() = $sseq =~ /-+/g) || 0;
				my $pid     =  (($iden/$ali_len) * 100);				
				my $qcov    =  (($ali_len/(length $qseq)) * 100); 				
				my $gaps    = $matches->{$sid}->{$from}->{GAPS};	
						
				
				#If all inclusion critreria are met (%id, %Qcov, num gaps, runs of gaps)
				
				if ( ($pid >= $opt->id) && ($qcov >= $opt->cov) && ($dashes <= $opt->gaps) && ($runs <= 1))
				{
					my @snts = split ("", $sseq);
					my @qnts = split ("", $qseq);
					my @mod_seq;
				
					for my $i (0..$#snts)
					{
						if ($snts[$i] =~ /\-/)
						{
							push @mod_seq, $qnts[$i]; 
						}
						else
						{
							push @mod_seq, $snts[$i]; 
						}
					}

					my $mod = join ("", @mod_seq);
					my $mod_aa = &gjoseqlib::translate_seq( $mod );
			
					my ($len, $strand);
					if ($from < $to)
					{
						$strand = "+";
						$len = ($to - $from) + 1;
					}
					elsif ($from > $to){
						$strand = "-";
						$len = ($from - $to) + 1;	
					}
			
					my $feature = {
						type        => $ft,
						contig      => $sid,
						aa_sequence => $mod_aa,
						location    => ([[$sid, $from, $strand, $len]]),
						product     => $anno,
						symbol      => $symbol,
						pssm        => ([[$fam, $name, $anno, "LowVan Transcript Edited Feature"]]),
					};
					
					push(@{$features{$ft}}, $feature);
				}
			
			
				# if the inclusion criteria are NOT met, but there is still a decent HSP
				# from the blast, we add it as a partial CDS that goes uncorrected.  No
				# protein translation is given.
				else
				{
					my $feature_type = "partial_cds";
					my ($len, $strand);
					if ($from < $to)
					{
						$strand = "+";
						$len = ($to - $from) + 1;
					}
					elsif ($from > $to){
						$strand = "-";
						$len = ($from - $to) + 1;	
					}

					my $lc_anno = lcfirst($anno);
				
					my $feature = {
						type        => $feature_type,
						contig      => $sid,
						location    => ([[$sid, $from, $strand, $len]]),
						product     => "Uncorrected "."$lc_anno"." encoding region",
						pssm        => ([[$fam, $name, $anno, "LowVan Transcript Edited Feature"]]),  #I don't actually use this but i left it there.
					};
					push(@{$features{$feature_type}}, $feature);
				}
			}
		}
	}	
	
	# Push features into the GTO
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
						-annotator           => 'LowVan Transcript Edited Feature',
						-protein_translation => $data->{aa_sequence},
						-function            => $data->{product},
						-family_assignments  => $data->{pssm},
						};				
					if (defined $data->{symbol} && $data->{symbol} ne '') 
					{
						$p->{-alias_pairs} = [[gene => $data->{symbol}]];
					}
					$genome_in->add_feature($p);
				}
				
				# Call a partial cds and do not add the AA seq if its a distant match.
				# No family assignment is generated
				elsif ($type eq 'partial_cds')
				{
					my $p = {
						-id	                 => $genome_in->new_feature_id($type),
						-type 	             => $type,
						-location 	         => $data->{location},
						-analysis_event_id 	 => $event_id,
						-annotator           => 'LowVan Transcript Edited Feature',
						-function            => $data->{product},
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
		print STDERR "Transcript edited features curated for: $fam, but none were found.\n";
		chdir($base);
		$genome_in->destroy_to_file($opt->output);			
 	}
 }


 else 
 {
 	print STDERR "No proteins from transcript editing for: $fam\n";
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
				my $gaps      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{gaps};
				my $sid       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{description}->[$k]->{title};
				my $hit_from  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_from};
				my $bit       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{bit_score};
				my $hit_to    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_to};
				my $qseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{qseq};
				my $hseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hseq};
								
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
					}				
				}
			}
		}
	}
	return $matches;
}

###########################################################













