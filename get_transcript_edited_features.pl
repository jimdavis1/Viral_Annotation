#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Time::HiRes 'gettimeofday';
use GenomeTypeObject;
use Getopt::Long::Descriptive;
use P3DataAPI;
use JSON::XS;
use File::Slurp;
use File::Path 'remove_tree';
use IPC::Run qw(run);
use Cwd;
use gjoseqlib;


my ($help, $tmp); 
my($opt, $usage) = describe_options("%c %o",
				    ["input|i=s"       => "Input GTO"],
				    ["output|o=s"      => "Output GTO"],
				    ["cov|c=i"         => "Minimum percent query coverage (D = 95)", { default => 95 }],
				    ["id|p=i"          => "Minimum percent identity  (D = 95)", { default => 95 }],
				    ["gaps|g=i"        => "Maximum number of allowable gaps (D = 2)", { default => 2 }],
		            ["threads|a=i"     => "Threads for the BLASTN (D = 24))", { default => 24 }],
				    ["json|j=s"        => "Full path to the JSON opts file", {default => "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json"}],
				    ["dir|d=s"         => "Full path to the directory of aligned transcripts", {default => "/home/jjdavis/bin/Viral_Annotation/Transcript-Editing"}],
				    ["tmp|t=s"         => "Declare name for temp dir (D = randomly named in cwd)"], 
				    ["help|h"          => "Show this help message"]);


print($usage->text), exit 0 if $opt->help;
die($usage->text) if @ARGV != 0;
if ($opt->tmp){$tmp = $opt->tmp;}
else{$tmp .= sprintf("%x", rand 16) for 1..20;}
my $dir = $opt->dir;

my $genome_in = GenomeTypeObject->create_from_file($opt->input);
$genome_in or die "Error reading and parsing input";
$genome_in->{features}->[0] or die "No features in GTO\n"; 

my $base = getcwd;


# We read the GTO to get the family 
my %pssm_fam;
for my $i (0 .. $#{$genome_in->{features}}) 
{
	if (($genome_in->{features}->[$i]->{type} =~ /CDS/) && ($genome_in->{features}->[$i]->{family_assignments}->[0]->[3] =~ /annotate_by_viral_pssm/))
	{
		$pssm_fam{$genome_in->{features}->[$i]->{family_assignments}->[0]->[0]}++;
	}
}
die "More than one viral family of PSSMs in GTO\n" if scalar(keys %pssm_fam) > 1;
my $fam = (keys %pssm_fam)[0];
$fam or die "GTO has no annotations from annotate_by_viral_pssm tool\n"; 


# Next we read the JSON to see if there are any transcript edited features
my $json      = decode_json(read_file($opt->json));
$genome_in or die "Error reading json protein feature data";

my @to_analyze;
foreach (keys %{$json->{$fam}->{features}})
{
	my $prot = $_;
	if ($json->{$fam}->{features}->{$prot}->{special} eq "transcript_edit")
	{
		my $anno = $json->{$fam}->{features}->{$prot}->{anno};
		my $feature_type = $json->{$fam}->{features}->{$prot}->{feature_type};
		push @to_analyze, ([$prot, $anno, $feature_type]);
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
  	 print STDERR "makeblastdb failed with rc=$?. Stdout:\n";
	}
	
	# create the GTO analysis event.
	my $event = {
 	   tool_name => "get_transcript_edited_features",
  	   execution_time => scalar gettimeofday,
	};
	my $event_id = $genome_in->add_analysis_event($event);


	#cycle through the transcript edited features and search for them one at a time with blastn.
	for my $i (0..$#to_analyze)
	{
		my $name = $to_analyze[$i][0];
		my $anno = $to_analyze[$i][1]; 
		my $ft   = $to_analyze[$i][2];
		
		print STDERR "\tAnalyzing $name\n\n"; 
		my $query = "$dir/$fam/$name.fasta"; 
		
		
			my $make_db2 = run("makeblastdb -dbtype nucl -in $query >/dev/null");

			if (!$make_db2)
			{
  				 print STDERR "makeblastdb failed with rc=$?. Stdout:\n";
			}

			my @blast_parms = (
		      "-query",         $query,
		      "-db",            $contigs, #### fix this in the original program
		      "-evalue",        0.5,
		      "-reward",          2,
		      "-penalty",        -3,
		      "-word_size",      28,
		      "-outfmt",         15,
		      "-soft_masking",   "false",
		      "-dust",           "no",
		      "-perc_identity",  $opt->id,
              "-qcov_hsp_perc",  $opt->cov,
		      "-num_threads",    $opt->threads);

		#print STDERR Dumper(\@blast_parms);
		my $do_blast = run(["blastn", @blast_parms], ">", "$name.json", "2>", "$name.blastn.stderr.txt");
		
		open (IN, "<$name.json"), or warn "Cannot open JSON BLASTn output file $name.json\n";
		my $results = decode_json(scalar read_file(\*IN));	
		close IN;
		
		#Gather in the best match.
		my $matches = best_blastn_match_by_loc($results, $opt->id, $opt->cov, $opt->gaps);  
	
		# Create the transcript edited protein feature by reading the subject and filling the 
		# gaps using the hand-curated sequences from the query. 	
		foreach (keys %$matches)
		{
			my $sid = $_; 
			foreach (keys %{$matches->{$sid}})
			{
				my $from = $_;
				my $to   = $matches->{$sid}->{$from}->{TO};
				my $sseq = $matches->{$sid}->{$from}->{HSEQ};
				my $qseq = $matches->{$sid}->{$from}->{QSEQ};
				my $dashes = ($sseq =~ tr/-//); 
				my $runs = $sseq =~ /(?=--)/g;

				#There can only be one gap with <= $opt->gaps number of dashes
				if (($dashes <= $opt->gaps) && ($runs <= 1))
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

				
					#okay, now we need to move the feature into the GTO
				
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
						pssm        => ([[$fam, $name, $anno, "get_transcript_edited_features"]]),
						};
					push(@{$features{$ft}}, $feature);
				}
			}
		}
	}

    for my $type (keys %features)
    {
		my $feats = $features{$type};
		my $n = @$feats;
		my $id_type = $type;
	
		for my $feature (@$feats)
		{
			my $p = {
					-id	                 => $genome_in->new_feature_id($id_type),
					-type 	             => $type,
					-location 	         => $feature->{location},
					-analysis_event_id 	 => $event_id,
					-annotator           => 'get_transcript_edited_features',
					-protein_translation => $feature->{aa_sequence},
					-function            => $feature->{product},
					-family_assignments  => $feature->{pssm},
					};
			$genome_in->add_feature($p);
		}
    }

	chdir ($base);
	remove_tree ($tmp);
	$genome_in->destroy_to_file($opt->output);
}			


else 
{
	print STDERR "No proteins from transcript editing for $fam\n";
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
#   This reads the -db formatted json output not the -subject [fasta] formatted version
#   The returned JSONs are slightly different 
#
#   usage:
#   $hash = best_blastn_match_by_loc($blastn_json, $Min%ID, $Min%Qcov, $max#gaps);
#   Note that the input %identity and %qcov are actually percentages and not fractions.
# 
#   The returned hash reference is in the following format:
#     
#    SubjectID->{HitFrom}->{TO}->{HitTo} 
#             ->{HitFrom}->{HSEQ}->{SubjectSeq}
#             ->{HitFrom}->{QSEQ}->{QuerySeq}
#             ->{HitFrom}->{BIT}->{BitScore}
# 
#----------------------------------------------------------
sub best_blastn_match_by_loc
{
	my ($blast, $min_pid, $min_qcov, $max_gaps)  =  @_;
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
				my $gaps      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{description}->[$k]->{gaps};
				my $sid       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{description}->[$k]->{title};
				my $hit_from  = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_from};
				my $bit       = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{bit_score};
				my $hit_to    = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hit_to};
				my $qseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{qseq};
				my $hseq      = $blast->{BlastOutput2}->[$i]->{report}->{results}->{search}->{hits}->[$j]->{hsps}->[$k]->{hseq};
								
				my ($pid, $qcov);
				if ($ident && $ali_len && $qlen) # this ensures that we got search results.
				{
					$pid  =  (($ident/$ali_len) * 100);
					$qcov =  (($ali_len/$qlen)  * 100); 
					
					if (($pid >= $min_pid) && ($qcov >= $min_qcov) && ($gaps <= $max_gaps))
					{
						# if we have already seen a match in this contig:
						if (exists $matches->{$sid}) 
						{							
							foreach (keys %{$matches->{$sid}})
							{
								my $loc = $_; 
								
								#if the match is in roughly in the same location as seen before, and it has a better bitscore, 
								#then they are considered the same match.  We delete the orignal and keep the one with the better bit score.							
								if ((abs($hit_from - $loc) < $ali_len) && ($bit > $matches->{$sid}->{$loc}->{BIT}))
								{
									delete $matches->{$sid}->{$loc};
									$matches->{$sid}->{$hit_from}->{BIT}   = $bit;
									$matches->{$sid}->{$hit_from}->{TO}    = $hit_to;
									$matches->{$sid}->{$hit_from}->{QSEQ}  = $qseq;
									$matches->{$sid}->{$hit_from}->{HSEQ}  = $hseq;
								}
							}
						}
						# if we havent seen a match in this location before add it:
						else  
						{
							$matches->{$sid}->{$hit_from}->{BIT}   = $bit;
							$matches->{$sid}->{$hit_from}->{TO}    = $hit_to;
							$matches->{$sid}->{$hit_from}->{QSEQ}  = $qseq;
							$matches->{$sid}->{$hit_from}->{HSEQ}  = $hseq;
						}
					}				
				}
			}
		}
	}
	return $matches;
}

###########################################################

























