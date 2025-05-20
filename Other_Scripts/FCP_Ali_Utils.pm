package FCP_Ali_Utils;

use strict;
use warnings;
use Math::Round;
use Cwd;
use Exporter qw(import);
use gjoseqlib;
use FCP_PSSM_Utils qw(create_pssm);

our @EXPORT_OK = qw(process_cluster_alignments process_alignment);

# Main function to process all cluster alignments
sub process_cluster_alignments 
{
    my ($seq_data, $options, $base) = @_;
    my $tmp = $seq_data->{tmp};
    
    # Get cluster files
    opendir(my $dir_h, "clusters") or die "Cannot open clusters directory: $!\n";
    my @clusterF = grep {$_ !~ /^\./} readdir($dir_h);
    closedir($dir_h);

    # Create report file
    open(my $report_fh, ">Curation_Report") or die "Cannot open Curation_Report for writing: $!\n";
    print $report_fh "Alignment\tAli_Len\tFirst_AA\tOriginal Seqs\tFinal Seqs\tLow_Occ_Cols_Cut\tN_Term_Cut\tC_Term_Cut\tDash_Starts_Removed\tCut_Cols\tSeq_w_End_Dashes_Removed\n";

    # Process each cluster file
    foreach my $file (@clusterF) 
    {
        my $ali_file = $file;
        $ali_file =~ s/\.fasta//g;
        
        print STDERR "\n\nAligning $file\n";
        
        # Align the sequences with MAFFT
        open(my $in_fh, "mafft --thread 24 --quiet --reorder $base/$tmp/clusters/$file |") or die "Cannot run MAFFT: $!\n";
        my @ali = &gjoseqlib::read_fasta($in_fh);
        close($in_fh);

        # Print the original alignment
        open(my $out_fh, ">alis/$ali_file.fa") or die "Cannot open alis/$ali_file.fa for writing: $!\n";
        &gjoseqlib::print_alignment_as_fasta($out_fh, @ali);
        close($out_fh);

        # Process the alignment
        my $ali_data = process_alignment(\@ali, $options);
        
        # Write to report file
        write_alignment_report($report_fh, $ali_file, $ali_data);
        
        # Write corrected alignment
        open($out_fh, ">corrected_alis/$ali_file.fa") or die "Cannot open corrected_alis/$ali_file.fa for writing: $!\n";
        &gjoseqlib::print_alignment_as_fasta($out_fh, @{$ali_data->{final_ali}});
        close($out_fh);
        
        # Create PSSM if enough sequences
        if ($ali_data->{n_seqs_final} >= $options->{min_seqs}) 
        {
            create_pssm_for_alignment($ali_file, $ali_data, $options, $tmp);
        }
    }
    
    close($report_fh);
    return 1;
}

# Process a single alignment
sub process_alignment 
{
    my ($ali, $options) = @_;
    
    # Extract options
    my $frac_dash = $options->{frac_dash};
    my $f_insert = $options->{f_insert};
    my $frac_n_term = $options->{frac_n_term};
    my $frac_c_term = $options->{frac_c_term};
    my $start_dash = $options->{start_dash};
    my $end_dash = $options->{end_dash};
    my $end_frac_occ = $options->{end_frac_occ};
    
    # Calculate original sequence count
    my $n_seqs_orig = scalar @$ali;
    
    # Remove sequences with too many dashes
    my @ali2;
    for my $i (0..$#$ali) 
    {
        my $len = length($ali->[$i][2]);
        my $dashes = $ali->[$i][2] =~ tr/-//;
        if (($dashes/$len) <= $frac_dash) 
        {
            push @ali2, $ali->[$i];
        }
    }
    
    my $n_seqs_internal_dash = scalar @ali2;
    
    # Pack the alignment (remove empty columns)
    my @ali3 = &gjoseqlib::pack_alignment(@ali2);
    
    # Process column occupancy
    my $n_seqs = scalar @ali3;
    my $min_n = round($f_insert * $n_seqs);
    my %col_sum;  # number of non-dash characters
    my $aa_sum = {}; # aa count
    
    # Count characters in each column
    for my $i (0..$#ali3) 
    {
        my $sequence = uc $ali3[$i][2];
        my @bases = split("", $sequence);
        for my $j (0..$#bases) 
        {
            if ($bases[$j] =~ /\w/) 
            {
                $aa_sum->{$j}->{$bases[$j]}++;  # only counts letters not dashes
                $col_sum{$j}++;
            }
        }
    }
    
    # Add low occupancy columns to exclusion hash
    my %exclude;
    foreach my $col (keys %col_sum) 
    {
        my $f_occ = ($col_sum{$col}/$n_seqs);
        if ($f_occ < $f_insert) 
        {
            print STDERR "Excluding $col\t$f_occ\n";
            $exclude{$col} = 1;
        }
    }
    
    my $low_occ = keys %exclude;  # count of low occupancy columns removed
    my $ali_len = length($ali3[0][2]);
    
    # Trim from the N-terminus
    my $n_term_trimmed = 0;
    my $c_term_trimmed = 0;
    
    for my $k (0..$ali_len) 
    {
        my $aa_resR = $aa_sum->{$k};
        next unless $aa_resR; # Skip if no data for this position
        
        my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys(%$aa_resR);
        next unless @sorted_res; # Skip if no residues
        
        my $mc_res = $sorted_res[0];
        my $cons = ($aa_sum->{$k}->{$mc_res}/$n_seqs_internal_dash);
        
        if ($cons < $frac_n_term) 
        {
            $exclude{$k} = 1;
            $n_term_trimmed++;
            print STDERR "Trimmed $k $mc_res, $aa_sum->{$k}->{$mc_res}, $cons, from N-terminal\n";
        } 
        else 
        {
            last;
        }
    }
    
    # Trim from the C-terminus
    for (my $l = ($ali_len - 1); $l >= 0; $l--) 
    {
        my $aa_resR = $aa_sum->{$l};
        next unless $aa_resR; # Skip if no data for this position
        
        my @sorted_res = sort {$aa_resR->{$b} <=> $aa_resR->{$a}} keys(%$aa_resR);
        next unless @sorted_res; # Skip if no residues
        
        my $mc_res = $sorted_res[0];
        my $cons = ($aa_sum->{$l}->{$mc_res}/$n_seqs_internal_dash);
        
        if ($cons < $frac_c_term) 
        {
            $exclude{$l} = 1;
            $c_term_trimmed++;
            print STDERR "Trimmed $l, $mc_res, $aa_sum->{$l}->{$mc_res}, $cons, from C-terminal\n";
        } 
        else 
        {
            last;
        }
    }
    
    # Apply exclusions and remove sequences starting with a dash if needed
    my @ali4;
    my $dash_starts_removed = 0;
    
    for my $i (0..$#ali3) 
    {
        my $id = $ali3[$i][0];
        my $anno = $ali3[$i][1];
        my @bases = split("", $ali3[$i][2]);
        my $string = "";
        
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
        elsif ($string !~ /^-/) 
        {
            push @ali4, ([$id, $anno, $string]);
        } 
        else 
        {
            $dash_starts_removed++;
        }
    }
    
    # Pack the alignment again
    my @ali5 = &gjoseqlib::pack_alignment(@ali4);
    
    # Get rid of identical sequences
    my $unique = {};
    for my $i (0..$#ali5) 
    {
        my $seq = uc $ali5[$i][2];
        $unique->{$seq}->{ID} = $ali5[$i][0];
        $unique->{$seq}->{ANNO} = $ali5[$i][1];
    }
    
    my @ali6;
    foreach my $seq (keys %$unique) 
    {
        push @ali6, ([$unique->{$seq}->{ID}, $unique->{$seq}->{ANNO}, $seq]);
    }
    
    # Handle sequences with spurious C-terminal gaps
    my %occ_cnt;
    my $len = length($ali6[0][2]);
    my $start_from = ($len - $end_dash);
    my $continue_cut = 1;
    
    # First pass: populate %occ_cnt for relevant columns at the end
    for my $i (0..$#ali6) 
    {
        my @bases = split("", $ali6[$i][2]);
        for my $j ($start_from..$#bases) 
        {
            if ($bases[$j] =~ /\w/) 
            {
                $occ_cnt{$j}++;
            }
        }
    }
    
    # Check if any end columns have too many dashes (might indicate conservation)
    my $nseqs = scalar @ali6;
    foreach my $col (sort {$a <=> $b} keys %occ_cnt) 
    {
        if (($occ_cnt{$col}/$nseqs) < (1 - $end_frac_occ)) 
        {
            $continue_cut = 0;
        }
    }
    
    my @ali7;
    if ($continue_cut) 
    {
        # Keep sequences that do not end in a dash
        for my $i (0..$#ali6) 
        {
            my $end = substr($ali6[$i][2], -$end_dash);
            if ($end =~ /\w/) 
            {
                push @ali7, $ali6[$i];
            }
        }
    } 
    else 
    {
        @ali7 = @ali6;
    }
    
    my @ali8 = &gjoseqlib::pack_alignment(@ali7);
    my $end_dash_removed = (scalar @ali6) - (scalar @ali7);
    
    # Final stats
    my $n_seqs_final = scalar @ali8;
    my $start_char = substr($ali8[0][2], 0, 1) if @ali8;
    my $final_ali_len = @ali8 ? length($ali8[0][2]) : 0;
    
    print STDERR "Original_Seqs = $n_seqs_orig, Final_Seqs = $n_seqs_final\n";
    print STDERR "Low Occ Cols Cut = $low_occ, N-term Cut = $n_term_trimmed, C-term Cut = $c_term_trimmed\t";
    print STDERR "Dash Starts Removed = $dash_starts_removed, Dash End Seqs Removed = $end_dash_removed\n\n";
    
    # Return all the data
    return {
        n_seqs_orig => $n_seqs_orig,
        n_seqs_final => $n_seqs_final,
        ali_len => $final_ali_len,
        start_char => $start_char,
        low_occ => $low_occ,
        n_term_trimmed => $n_term_trimmed,
        c_term_trimmed => $c_term_trimmed,
        dash_starts_removed => $dash_starts_removed,
        end_dash_removed => $end_dash_removed,
        exclude => \%exclude,
        final_ali => \@ali8
    };
}

# Write alignment processing report
sub write_alignment_report 
{
    my ($report_fh, $ali_file, $ali_data) = @_;
    
    print $report_fh join("\t", 
        $ali_file,
        $ali_data->{ali_len},
        $ali_data->{start_char},
        $ali_data->{n_seqs_orig},
        $ali_data->{n_seqs_final},
        $ali_data->{low_occ},
        $ali_data->{n_term_trimmed},
        $ali_data->{c_term_trimmed},
        $ali_data->{dash_starts_removed},
        $ali_data->{end_dash_removed},
        join(",", sort {$a <=> $b} keys %{$ali_data->{exclude}})
    ), "\n";
    
    return 1;
}

# Create PSSM for an alignment
sub create_pssm_for_alignment 
{
    my ($ali_file, $ali_data, $options, $tmp) = @_;
    
    my $pssm_file = $ali_file;
    if ($options->{pssm_prefix}) 
    {
        $pssm_file = "$options->{pssm_prefix}.$ali_file";
    }
    
    # Set up PSSM options
    my $pssm_options = {
        outfile => "pssms/$tmp.pssm",
        final_file => "pssms/$pssm_file.pssm",
        process_titles => 1,
        cleanup => 1,
        tmp => $tmp,
        pssm_prefix => $options->{pssm_prefix}
    };
    
    # Create the PSSM using PSSMUtils
    create_pssm($ali_data->{final_ali}, $pssm_options);
    
    return 1;
}

1; # End of module