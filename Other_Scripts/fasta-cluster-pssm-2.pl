#!/usr/bin/env perl
use strict;
use warnings;
use lib '.';  # Add current directory to module search path
use Getopt::Long;
use Cwd;
use Data::Dumper;
use Math::Round;

# Import our modules with underscores instead of hyphens
use FCP_Main_Utils qw(parse_command_line process_sequences run_mmseqs);
use FCP_Ali_Utils qw(process_cluster_alignments);
use FCP_PSSM_Utils qw(create_pssm);
use FCP_Nterm_Utils qw(evaluate_nterm_conservation recluster_by_nterm);

# Import required external modules
use gjoseqlib;
use BlastInterface;

# Main script execution
print STDERR "Parsing command line...\n";
my $options = parse_command_line();
my $tmp_dir = $options->{tmp};

# Generate random tmp directory name if not specified
if (!$tmp_dir) 
{
    $tmp_dir = "";
    $tmp_dir .= sprintf("%x", rand 16) for 1..10;
    $options->{tmp} = $tmp_dir; # Update options hash with the new tmp_dir
}
print STDERR "Using temp directory: $tmp_dir\n";

# Read proteins from stdin
print STDERR "Reading protein sequences from stdin...\n";
my @prots = &gjoseqlib::read_fasta();
print STDERR "Read ", scalar(@prots), " sequences\n";

# Create the tempdir and go there
my $base = getcwd;
print STDERR "Current directory: $base\n";
mkdir($tmp_dir) or die "Cannot create directory $tmp_dir: $!\n";
chdir($tmp_dir) or die "Cannot change to directory $tmp_dir: $!\n";
print STDERR "Changed to directory: ", getcwd, "\n";

# Process sequences and create unique set
print STDERR "Processing sequences...\n";
my $seq_data = process_sequences(\@prots, $options);
print STDERR "Sequence processing complete\n";

# Run MMSeqs clustering
print STDERR "Running MMSeqs clustering...\n";
run_mmseqs($seq_data, $options);
print STDERR "MMSeqs clustering complete\n";

# Process cluster alignments
print STDERR "Processing cluster alignments...\n";
process_cluster_alignments($seq_data, $options, $base);
print STDERR "Alignment processing complete\n";

# Create directories needed for N-terminal evaluation
mkdir("alis-for-reclustering") unless -d "alis-for-reclustering";
mkdir("reclustering") unless -d "reclustering";

# Evaluate N-terminal conservation and recluster if needed
# This is now enabled by default
if ($options->{evaluate_nterm}) 
{
    print STDERR "Evaluating N-terminal conservation...\n";
    my $results = evaluate_nterm_conservation($options);
    
    if ($results->{alignments_to_recluster} && @{$results->{alignments_to_recluster}}) 
    {
        print STDERR "Found ", scalar(@{$results->{alignments_to_recluster}}), " alignments with poor N-terminal conservation\n";
        print STDERR "Reclustering these alignments with identity threshold: $options->{n_term_cluster_id}...\n";
        recluster_by_nterm($results->{alignments_to_recluster}, $options, $base);
    } 
    else 
    {
        print STDERR "All alignments have good N-terminal conservation\n";
    }
}

print STDERR "All processing complete\n";
exit(0);
