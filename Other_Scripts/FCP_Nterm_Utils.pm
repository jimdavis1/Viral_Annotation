package FCP_Nterm_Utils;

use strict;
use warnings;
use Exporter qw(import);
use File::Basename;
use File::Copy;
use Cwd qw(abs_path);
use Data::Dumper;
use gjoseqlib;
use FCP_PSSM_Utils qw(create_pssm);

our @EXPORT_OK = qw(evaluate_nterm_conservation recluster_by_nterm);

# Constants
use constant NTERM_LENGTH => 25; # Number of amino acids to analyze from N-terminus

# Evaluate N-terminal conservation in all corrected alignments
sub evaluate_nterm_conservation 
{
    my ($options) = @_;
    my $p_nterm = $options->{p_nterm} / 100.0; # Convert percentage to fraction
    my $n_nterm = $options->{n_nterm};
    
    print STDERR "Evaluating N-terminal conservation (threshold: $options->{p_nterm}%, min columns: $n_nterm)...\n";
    
    # Get list of all alignments in corrected_alis directory
    opendir(my $dh, "corrected_alis") or die "Cannot open corrected_alis directory: $!\n";
    my @alignment_files = grep { /\.fa$/ } readdir($dh);
    closedir($dh);
    
    print STDERR "Found ", scalar(@alignment_files), " alignments to evaluate\n";
    
    # Create alis-for-reclustering directory if it doesn't exist
    mkdir("alis-for-reclustering") unless -d "alis-for-reclustering";
    
    # Open log file
    open(my $log_fh, ">nterm_evaluation.log") or die "Cannot open nterm_evaluation.log for writing: $!\n";
    print $log_fh "Alignment\tTotal_Columns\tNterm_Columns\tPoor_Columns\tStatus\n";
    
    my @alignments_to_recluster;
    
    # Process each alignment
    foreach my $ali_file (@alignment_files) 
    {
        print STDERR "Evaluating N-terminal conservation for $ali_file...\n";
        
        # Read the alignment
        open(my $ali_fh, "<corrected_alis/$ali_file") or die "Cannot open corrected_alis/$ali_file: $!\n";
        my @ali = &gjoseqlib::read_fasta($ali_fh);
        close($ali_fh);
        
        # Skip if alignment is empty
        if (!@ali) 
        {
            print $log_fh "$ali_file\t0\t0\t0\tSkipped (empty)\n";
            next;
        }
        
        my $n_seqs = scalar @ali;
        
        # Extract N-terminal regions (first NTERM_LENGTH amino acids)
        my @nterm_ali;
        my $max_nterm_length = 0;
        
        foreach my $seq (@ali) 
        {
            my $seq_data = $seq->[2];
            my $nterm_seq = substr($seq_data, 0, NTERM_LENGTH);
            push @nterm_ali, [$seq->[0], $seq->[1], $nterm_seq];
            
            my $len = length($nterm_seq);
            $max_nterm_length = $len if $len > $max_nterm_length;
        }
        
        # Skip if N-terminal region is too short
        if ($max_nterm_length < 5) 
        {
            print $log_fh "$ali_file\t", length($ali[0][2]), "\t$max_nterm_length\t0\tSkipped (N-terminal too short)\n";
            next;
        }
        
        # Count amino acid frequencies for each column
        my %col_stats;
        
        for my $i (0..$#nterm_ali) 
        {
            my $seq = $nterm_ali[$i][2];
            my @residues = split('', $seq);
            
            for my $j (0..$#residues) 
            {
                next if $residues[$j] eq '-'; # Skip gaps
                $col_stats{$j}{$residues[$j]}++;
            }
        }
        
        # Identify poorly conserved columns
        my @poor_columns;
        
        for my $j (0..$max_nterm_length-1) 
        {
            next unless exists $col_stats{$j};
            
            # Count total non-gap residues in this column
            my $total = 0;
            foreach my $aa (keys %{$col_stats{$j}}) 
            {
                $total += $col_stats{$j}{$aa};
            }
            
            next if $total < 3; # Skip columns with too few residues
            
            # Find most frequent amino acid and its fraction
            my $max_aa = '';
            my $max_count = 0;
            
            foreach my $aa (keys %{$col_stats{$j}}) 
            {
                if ($col_stats{$j}{$aa} > $max_count) 
                {
                    $max_count = $col_stats{$j}{$aa};
                    $max_aa = $aa;
                }
            }
            
            my $conservation = $max_count / $total;
            
            # If conservation below threshold, mark as a poor column
            if ($conservation < $p_nterm) 
            {
                push @poor_columns, $j;
            }
        }
        
        # Determine if we need to recluster based on the number of poor columns
        my $status = "OK";
        if (scalar(@poor_columns) >= $n_nterm) 
        {
            $status = "Needs reclustering";
            push @alignments_to_recluster, {
                file => $ali_file,
                alignment => \@ali,
                nterm_alignment => \@nterm_ali,
                poor_columns => \@poor_columns
            };
            
            # Move the original alignment to alis-for-reclustering directory (not copy)
            rename("corrected_alis/$ali_file", "alis-for-reclustering/$ali_file") or 
                die "Cannot move corrected_alis/$ali_file to alis-for-reclustering/$ali_file: $!\n";
        }
        
        # Log result
        print $log_fh "$ali_file\t", length($ali[0][2]), "\t$max_nterm_length\t", scalar(@poor_columns), "\t$status\n";
        
        if ($status eq "Needs reclustering") 
        {
            print STDERR "  $ali_file: Found ", scalar(@poor_columns), " poorly conserved N-terminal columns, marking for reclustering\n";
        } 
        else 
        {
            print STDERR "  $ali_file: N-terminal conservation is acceptable\n";
        }
    }
    
    close($log_fh);
    
    print STDERR "N-terminal evaluation complete, ", scalar(@alignments_to_recluster), " alignments need reclustering\n";
    
    return {
        alignments_to_recluster => \@alignments_to_recluster
    };
}

# Recluster alignments based on N-terminal conservation
sub recluster_by_nterm 
{
    my ($alignments, $options, $base_dir) = @_;
    
    print STDERR "Reclustering ", scalar(@$alignments), " alignments...\n";
    
    # Create directories if they don't exist
    mkdir("reclustered_alis") unless -d "reclustered_alis";
    mkdir("reclustering") unless -d "reclustering";
    
    # Open log file for reclustering
    open(my $log_fh, ">nterm_reclustering.log") or die "Cannot open nterm_reclustering.log for writing: $!\n";
    print $log_fh "Original_Alignment\tClusters\tNew_Alignments\tLeftover_Seqs\n";
    
    # Process each alignment needing reclustering
    foreach my $ali_data (@$alignments) 
    {
        my $ali_file = $ali_data->{file};
        my $ali_base = basename($ali_file, '.fa');
        
        print STDERR "Reclustering $ali_file...\n";
        
        # Remove any existing PSSM for this alignment
        my $old_pssm_file = $ali_base;
        if ($options->{pssm_prefix}) {
            $old_pssm_file = "$options->{pssm_prefix}.$ali_base";
        }
        
        my $old_pssm_path = "pssms/$old_pssm_file.pssm";
        if (-f $old_pssm_path) {
            print STDERR "  Removing old PSSM: $old_pssm_path\n";
            unlink($old_pssm_path) or 
                print STDERR "  Warning: Could not remove old PSSM file $old_pssm_path: $!\n";
        }
        
        # Create temp directory for this alignment inside reclustering
        my $nterm_dir = "reclustering/nterm_$ali_base";
        mkdir($nterm_dir) or die "Cannot create directory $nterm_dir: $!\n";
        
        # Extract N-terminal regions for clustering
        open(my $nterm_fh, ">$nterm_dir/nterm.fasta") or die "Cannot open $nterm_dir/nterm.fasta for writing: $!\n";
        foreach my $seq (@{$ali_data->{nterm_alignment}}) 
        {
            &gjoseqlib::print_alignment_as_fasta($nterm_fh, $seq);
        }
        close($nterm_fh);
        
        # Run MMSeqs on the N-terminal alignment with custom identity threshold
        print STDERR "  Running MMSeqs on N-terminal region...\n";
        my $mmseqs_prefix = "$nterm_dir/nterm";
        
        # Use the n_term_cluster_id parameter with a default of 0.7
        my $n_term_cluster_id = defined($options->{n_term_cluster_id}) ? $options->{n_term_cluster_id} : 0.7;
        my $mmc = $options->{mmc};
        my $cmd = "mmseqs easy-cluster $nterm_dir/nterm.fasta $mmseqs_prefix tmp --min-seq-id $n_term_cluster_id -c $mmc --cov-mode 0 >/dev/null";
        print STDERR "  Command: $cmd\n";
        
        my $ret = system $cmd;
        if ($ret != 0) 
        {
            print STDERR "  Warning: MMSeqs execution failed with return code $ret: $!\n";
            print $log_fh "$ali_file\t0\t0\t0\n";
            next;
        }
        
        # Read the all_seqs fasta file
        print STDERR "  Reading $mmseqs_prefix\_all_seqs.fasta...\n";
        open(my $in_fh, "<$mmseqs_prefix\_all_seqs.fasta") or die "Cannot open $mmseqs_prefix\_all_seqs.fasta: $!\n";
        my $seqH = {};
        my @seqs = &gjoseqlib::read_fasta($in_fh);
        close($in_fh);
        
        foreach my $seq (@seqs) 
        {
            $seqH->{$seq->[0]} = {
                ANNO => $seq->[1],
                SEQ => $seq->[2]
            };
        }
        
        # Read the cluster TSV file
        print STDERR "  Reading $mmseqs_prefix\_cluster.tsv...\n";
        open($in_fh, "<$mmseqs_prefix\_cluster.tsv") or die "Cannot open $mmseqs_prefix\_cluster.tsv: $!\n";
        my %clusters;
        while (<$in_fh>) 
        {
            chomp;
            my ($id, $memb) = split /\t/;
            push @{$clusters{$id}}, $memb;
        }
        close($in_fh);
        
        print STDERR "  Found ", scalar(keys %clusters), " clusters\n";
        
        # Create a mapping from sequence ID to original full sequence
        my %id_to_orig;
        foreach my $seq (@{$ali_data->{alignment}}) 
        {
            $id_to_orig{$seq->[0]} = $seq;
        }
        
        # Process clusters and create new alignments
        my $min_seqs = $options->{min_seqs};
        my $new_alignments = 0;
        my $leftover_seqs = 0;
        
        # Prepare leftover sequences file
        open(my $leftover_fh, ">$nterm_dir/leftover_seqs.fa") or die "Cannot open $nterm_dir/leftover_seqs.fa for writing: $!\n";
        
        # Process each cluster
        my $cluster_num = 1;
        my @new_alignment_files; # Store new alignment filenames for PSSM creation
        
        foreach my $cluster_id (sort { scalar(@{$clusters{$b}}) <=> scalar(@{$clusters{$a}}) } keys %clusters) 
        {
            my @members = @{$clusters{$cluster_id}};
            
            if (scalar(@members) >= $min_seqs) 
            {
                # Create a new alignment file
                my $out_file = "reclustered_alis/${ali_base}_${cluster_num}.fa";
                print STDERR "  Creating new alignment $out_file with ", scalar(@members), " sequences\n";
                
                open(my $out_fh, ">$out_file") or die "Cannot open $out_file for writing: $!\n";
                
                # Write full sequences to the new alignment file
                foreach my $seq_id (@members) 
                {
                    if (exists $id_to_orig{$seq_id}) 
                    {
                        &gjoseqlib::print_alignment_as_fasta($out_fh, $id_to_orig{$seq_id});
                    } 
                    else 
                    {
                        print STDERR "  Warning: Sequence ID $seq_id not found in original alignment\n";
                    }
                }
                
                close($out_fh);
                push @new_alignment_files, {
                    file => $out_file,
                    base => "${ali_base}_${cluster_num}"
                };
                $cluster_num++;
                $new_alignments++;
                
                # Align the new file with MAFFT
                my $mafft_cmd = "mafft --thread 24 --quiet --reorder $out_file > ${out_file}.aligned";
                print STDERR "  Running MAFFT on $out_file\n";
                my $mafft_ret = system $mafft_cmd;
                
                if ($mafft_ret == 0) 
                {
                    # Replace the original file with the aligned version
                    rename("${out_file}.aligned", $out_file) or
                        die "Cannot rename ${out_file}.aligned to $out_file: $!\n";
                } 
                else 
                {
                    print STDERR "  Warning: MAFFT execution failed for $out_file\n";
                }
            } 
            else 
            {
                # Add to leftover sequences
                foreach my $seq_id (@members) 
                {
                    if (exists $id_to_orig{$seq_id}) 
                    {
                        &gjoseqlib::print_alignment_as_fasta($leftover_fh, $id_to_orig{$seq_id});
                        $leftover_seqs++;
                    }
                }
            }
        }
        
        close($leftover_fh);
        
        # If there are leftover sequences, add them to the global leftover file
        if ($leftover_seqs > 0) 
        {
            print STDERR "  Adding $leftover_seqs leftover sequences to Leftover_Seqs.aa\n";
            
            # Read the leftover sequences
            open(my $leftovers, "<$nterm_dir/leftover_seqs.fa") or die "Cannot open $nterm_dir/leftover_seqs.fa: $!\n";
            my @leftover_data = &gjoseqlib::read_fasta($leftovers);
            close($leftovers);
            
            # Append to the global leftover file
            open(my $global_leftover, ">>Leftover_Seqs.aa") or die "Cannot open Leftover_Seqs.aa for appending: $!\n";
            foreach my $seq (@leftover_data) 
            {
                &gjoseqlib::print_alignment_as_fasta($global_leftover, $seq);
            }
            close($global_leftover);
        }
            
        # Log results
        print $log_fh "$ali_file\t", scalar(keys %clusters), "\t$new_alignments\t$leftover_seqs\n";
        
        # Generate PSSM for each new alignment
        foreach my $new_ali (@new_alignment_files) 
        {
            print STDERR "  Creating PSSM for $new_ali->{file}...\n";
            generate_pssm_for_alignment($new_ali->{file}, $new_ali->{base}, $options);
        }
    }
    
    close($log_fh);
    
    # Re-align leftover sequences with MAFFT
    if (-s "Leftover_Seqs.aa") 
    {
        print STDERR "Re-aligning leftover sequences...\n";
        my $ret = system "mafft --thread 24 --quiet --reorder Leftover_Seqs.aa > Leftover_Seqs.fa";
        
        if ($ret != 0) 
        {
            print STDERR "Warning: MAFFT execution failed for Leftover_Seqs.aa: $!\n";
        }
    }
    
    print STDERR "Reclustering complete\n";
    return 1;
}

# Generate PSSM for an alignment
sub generate_pssm_for_alignment 
{
    my ($alignment_file, $base_name, $options) = @_;
    
    # Read the alignment
    open(my $ali_fh, "<$alignment_file") or die "Cannot open $alignment_file for reading: $!\n";
    my @ali = &gjoseqlib::read_fasta($ali_fh);
    close($ali_fh);
    
    # Skip if alignment is empty
    if (!@ali) {
        print STDERR "  Warning: Alignment in $alignment_file is empty, skipping PSSM generation\n";
        return;
    }
    
    # Set output file name
    my $pssm_file = $base_name;
    if ($options->{pssm_prefix}) {
        $pssm_file = "$options->{pssm_prefix}.$base_name";
    }
    
    # Set up PSSM options
    my $pssm_options = {
        outfile => "pssms/tmp_$base_name.pssm",
        final_file => "pssms/$pssm_file.pssm",
        process_titles => 1,
        cleanup => 1,
        tmp => $options->{tmp},
        pssm_prefix => $options->{pssm_prefix}
    };
    
    # Create the PSSM
    print STDERR "  Generating PSSM: pssms/$pssm_file.pssm\n";
    create_pssm(\@ali, $pssm_options);
    
    return 1;
}

1; # End of module