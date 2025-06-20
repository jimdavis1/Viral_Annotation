package FCP_Main_Utils;

use strict;
use warnings;
use Getopt::Long;
use Exporter qw(import);
use gjoseqlib;
use Data::Dumper;

our @EXPORT_OK = qw(parse_command_line process_sequences run_mmseqs);

# Parse command line arguments and return options hash
sub parse_command_line 
{
    my $usage = 'fasta-cluster-pssm.pl [parms] <Protein-Sequences.fasta
			-h help
			-a "annotation string" (this also becomes the title in the pssm). 
			-m min number of sequences for building pssm [d = 5]
			-f delete columns with less than frac bases (d = 0.33)
			-n delete columns from the n-terminus until >= frac identity (d = 0.75)
			-c delete columns from the c-terminus until >= frac identity (d = 0)
			-x include non-standard amino acids, by default it removes the 
			   entire sequence if it has ambiguous (B|J|X|Z) characters.
			-i include identical sequences [d = keep only unique]
			-d keep sequences beginning with a dash [d = throw away]
			-e cutoff for the number of dashes at the end of a sequence to be seen 
			   before throwing away sequence with a gap at the end [d = 3]
			-efo end fraction occupancy, this goes with -e.  This dictates the fraction of dashes in the -e 
			     end columns in order to allow cutting [d = 0.15] (if the fraction is too high, you might 
			     be cutting something that is conserved)
			-fd maximum fraction of dashes allowable for keeping a sequence (d = 0.20)
			-mi mmseqs identity [d = 0.8]
			-mc mmseqs coverage [d = 0.8]
			-k keep temp fasta (pre-aligned)
			-p "pssm prefix" (default = nothing).
			-tmp name temp fasta file 
			-p_nterm percent identity minimum for reclustering the columns in the first 25aas (default = 65%)
			-n_nterm number of bad columns required for the n-terminal reclustering procedure (default = 3)
			-nterm_eval_length number of N-terminal columns to check for poor conservation (default = 10)
			-nterm-mmseq-id identity threshold for N-terminal reclustering (default = 0.7)
			-no_nterm_eval disable N-terminal evaluation and reclustering
';

    my ($help, $keep_x, $keep_i, $keep_tmp, $start_dash, $annotation, $pssm_prefix, $no_nterm_eval);

    my $end_frac_occ = 0.15;
    my $end_dash = 3;
    my $min_seqs = 5;
    my $f_insert = 0.33;
    my $frac_n_term = 0.75;
    my $frac_c_term = 0;
    my $frac_dash = 0.2;
    my $mmi = 0.8;
    my $mmc = 0.8;
    my $tmp;
    my $p_nterm = 65; # Percentage for N-terminal conservation
    my $n_nterm = 3;  # Number of columns required
    my $nterm_eval_length = 10; # Default length for N-terminal evaluation
    my $n_term_cluster_id = 0.7; # Default for N-terminal reclustering

    my $opts = GetOptions(
        'h'     => \$help,
        'a=s'   => \$annotation,
        'm=s'   => \$min_seqs,
        'f=f'   => \$f_insert,
        'e=i'   => \$end_dash,
        'efo=f' => \$end_frac_occ,
        'n=f'   => \$frac_n_term,
        'c=f'   => \$frac_c_term,
        'fd=f'  => \$frac_dash,
        'x'     => \$keep_x,
        'i'     => \$keep_i,
        'k'     => \$keep_tmp,
        'mi=f'  => \$mmi,
        'mc=f'  => \$mmc,
        'tmp=s' => \$tmp,
        'p=s'   => \$pssm_prefix,
        'd'     => \$start_dash,
        'p_nterm=i' => \$p_nterm,
        'n_nterm=i' => \$n_nterm,
        'nterm_eval_length=i' => \$nterm_eval_length,
        'nterm-mmseq-id=f' => \$n_term_cluster_id, # New parameter for N-terminal reclustering identity
        'no_nterm_eval' => \$no_nterm_eval
    );

    if ($help) 
    {
        die "$usage\n";
    }
    
    print STDERR "Command line options parsed:\n";
    print STDERR "  annotation: ", ($annotation || "not specified"), "\n";
    print STDERR "  tmp: ", ($tmp || "not specified"), "\n";
    print STDERR "  pssm_prefix: ", ($pssm_prefix || "not specified"), "\n";
    
    # N-terminal evaluation is enabled by default, only disabled if no_nterm_eval is specified
    my $evaluate_nterm = $no_nterm_eval ? 0 : 1;
    
    if ($evaluate_nterm)
    {
        print STDERR "  N-terminal evaluation is enabled (p_nterm: $p_nterm%, n_nterm: $n_nterm columns)\n";
        print STDERR "  N-terminal evaluation length: $nterm_eval_length amino acids\n";
        print STDERR "  N-terminal reclustering identity threshold: $n_term_cluster_id\n";
    }
    else
    {
        print STDERR "  N-terminal evaluation is disabled\n";
    }

    return {
        keep_x        => $keep_x,
        keep_i        => $keep_i,
        keep_tmp      => $keep_tmp,
        start_dash    => $start_dash,
        annotation    => $annotation,
        pssm_prefix   => $pssm_prefix,
        end_frac_occ  => $end_frac_occ,
        end_dash      => $end_dash,
        min_seqs      => $min_seqs,
        f_insert      => $f_insert,
        frac_n_term   => $frac_n_term,
        frac_c_term   => $frac_c_term,
        frac_dash     => $frac_dash,
        mmi           => $mmi,
        mmc           => $mmc,
        tmp           => $tmp,
        p_nterm       => $p_nterm,
        n_nterm       => $n_nterm,
        nterm_eval_length => $nterm_eval_length,
        n_term_cluster_id => $n_term_cluster_id,
        evaluate_nterm => $evaluate_nterm
    };
}

# Process input sequences and create temporary fasta file
sub process_sequences 
{
    my ($prots, $options) = @_;
    my $tmp = $options->{tmp};
    
    print STDERR "Processing sequences with options:\n", Dumper($options), "\n";
    print STDERR "Input sequences: ", scalar(@$prots), "\n";
    
    # Get options from the options hash
    my $keep_i = $options->{keep_i} || 0;  # Initialize to 0 if undefined
    my $keep_x = $options->{keep_x} || 0;  # Initialize to 0 if undefined
    my $annotation = $options->{annotation};
    
    # Prepare for sequence processing
    my @prots2;
    my $uniq = {};
    
    print STDERR "Preparing sequences (keep_i=$keep_i, keep_x=$keep_x)...\n";
    
    if ($keep_i) 
    {
        # For very large datasets, use a more efficient approach to sorting
        print STDERR "Keeping identical sequences. Sorting by length...\n";
        # Calculate lengths once to avoid repeated calculation
        my @seq_with_len;
        foreach my $prot (@$prots) 
        {
            push @seq_with_len, [$prot->[0], $prot->[1], $prot->[2], length($prot->[2])];
        }
        
        print STDERR "Sorting ", scalar(@seq_with_len), " sequences...\n";
        
        # Sort by the pre-calculated length
        @prots2 = map { [$_->[0], $_->[1], $_->[2]] } 
                  sort { $b->[3] <=> $a->[3] } 
                  @seq_with_len;
        
        print STDERR "Sorted ", scalar(@prots2), " sequences by length.\n";
    } 
    else 
    {
        if ($keep_x) 
        {
            print STDERR "Including non-standard amino acids.\n";
            foreach my $prot (@$prots) 
            {
                $uniq->{$prot->[2]}->{ID}   = $prot->[0];
                $uniq->{$prot->[2]}->{ANNO} = $prot->[1];
            }
        } 
        else 
        {
            print STDERR "Removing sequences with ambiguous residues.\n";
            foreach my $prot (@$prots) 
            {
                unless ($prot->[2] =~ /(B|J|X|Z)/i) 
                {
                    $uniq->{$prot->[2]}->{ID}   = $prot->[0];
                    $uniq->{$prot->[2]}->{ANNO} = $prot->[1];
                }
            }
        }

        print STDERR "Found ", scalar(keys %$uniq), " unique sequences.\n";
        print STDERR "Building sequence array...\n";
        
        # Build array with length for sorting
        my @seq_with_len;
        foreach my $seq (keys %$uniq) 
        {
            push @seq_with_len, [$uniq->{$seq}->{ID}, $uniq->{$seq}->{ANNO}, $seq, length($seq)];
        }
        
        print STDERR "Sorting ", scalar(@seq_with_len), " unique sequences...\n";
        
        # Sort by the pre-calculated length
        @prots2 = map { [$_->[0], $_->[1], $_->[2]] } 
                  sort { $b->[3] <=> $a->[3] } 
                  @seq_with_len;
        
        print STDERR "Sorted ", scalar(@prots2), " unique sequences by length.\n";
    }

    # Print fasta to temp
    print STDERR "Writing sequences to $tmp.fasta...\n";
    open(my $out_fh, ">$tmp.fasta") or die "Cannot open $tmp.fasta for writing: $!\n";

    for my $i (0..$#prots2) 
    {
        my $id = $prots2[$i][0];
        $id =~ s/\|[^|]*$// if $id =~ tr/|// > 1;
        if ($annotation) 
        {
            &gjoseqlib::print_alignment_as_fasta($out_fh, ([$id, $annotation, $prots2[$i][2]]));
        } 
        else 
        {
            &gjoseqlib::print_alignment_as_fasta($out_fh, ([$id, $prots2[$i][1], $prots2[$i][2]]));
        }
    }
    
    close($out_fh);
    print STDERR "Wrote ", scalar(@prots2), " sequences to $tmp.fasta\n";
    
    return {
        prots => \@prots2,
        tmp => $tmp
    };
}

# Run MMSeqs easy-cluster on the sequence data
sub run_mmseqs 
{
    my ($seq_data, $options) = @_;
    my $tmp = $seq_data->{tmp};
    
    my $mmi = $options->{mmi};
    my $mmc = $options->{mmc};
    
    print STDERR "Running mmseqs with parameters: min-seq-id=$mmi, coverage=$mmc\n";
    
    # Run mmseqs easy-cluster
    print STDERR "Command: mmseqs easy-cluster $tmp.fasta $tmp.mmseq tmp --min-seq-id $mmi -c $mmc --cov-mode 0\n";
    my $ret = system "mmseqs easy-cluster $tmp.fasta $tmp.mmseq tmp --min-seq-id $mmi -c $mmc --cov-mode 0 >/dev/null";
    
    if ($ret != 0) 
    {
        die "MMSeqs execution failed with return code $ret: $!\n";
    }
    
    print STDERR "MMSeqs clustering completed successfully\n";
    
    # Process the clusters
    process_mmseqs_clusters($seq_data, $options);
    
    return 1;
}

# Process MMSeqs clustering results
sub process_mmseqs_clusters 
{
    my ($seq_data, $options) = @_;
    my $tmp = $seq_data->{tmp};
    my $min_seqs = $options->{min_seqs};
    
    print STDERR "Processing MMSeqs clusters (min_seqs=$min_seqs)...\n";
    
    # Read the all_seqs fasta file
    print STDERR "Reading $tmp.mmseq_all_seqs.fasta...\n";
    open(my $in_fh, "<$tmp.mmseq_all_seqs.fasta") or die "Cannot open all_seqs fasta file from mmseqs: $!\n";
    my $seqH = {};
    my @seqs = &gjoseqlib::read_fasta($in_fh);
    close($in_fh);
    
    print STDERR "Read ", scalar(@seqs), " sequences from $tmp.mmseq_all_seqs.fasta\n";
    
    foreach my $i (0..$#seqs) 
    {
        $seqH->{$seqs[$i][0]}->{ANNO} = $seqs[$i][1];
        $seqH->{$seqs[$i][0]}->{SEQ}  = $seqs[$i][2];
    }

    # Process the cluster TSV file
    print STDERR "Reading $tmp.mmseq_cluster.tsv...\n";
    open($in_fh, "<$tmp.mmseq_cluster.tsv") or die "Couldn't open tsv file from mmseqs: $!\n";
    my %clusters;
    while (<$in_fh>) 
    {
        chomp;
        my ($id, $memb) = split /\t/;
        push @{$clusters{$id}}, $memb;
    }
    close($in_fh);
    
    print STDERR "Found ", scalar(keys %clusters), " clusters\n";

    # Create cluster files
    # First make sure the clusters directory exists
    mkdir("clusters") unless -d "clusters";
    
    my $cluster_count = 1;
    open(my $out2_fh, ">Leftover_Seqs.aa") or die "Cannot open Leftover_Seqs.aa for writing: $!\n";

    my $large_clusters = 0;
    my $small_clusters = 0;
    my $leftover_seqs = 0;
    
    foreach my $cluster_id (sort { scalar(@{$clusters{$b}}) <=> scalar(@{$clusters{$a}}) } keys %clusters) 
    {
        my @array = @{$clusters{$cluster_id}};
        
        if ((scalar @array) > $min_seqs) 
        {
            print STDERR "Creating cluster file clusters/$cluster_count.fasta with ", scalar(@array), " sequences\n";
            open(my $out_fh, ">clusters/$cluster_count.fasta") or die "Cannot open clusters/$cluster_count.fasta for writing: $!\n";
            foreach my $seq_id (@array) 
            {
                unless (exists $seqH->{$seq_id}) 
                {
                    print STDERR "Warning: Sequence ID $seq_id not found in mmseq_all_seqs.fasta\n";
                    next;
                }
                &gjoseqlib::print_alignment_as_fasta($out_fh, ([$seq_id, $seqH->{$seq_id}->{ANNO}, $seqH->{$seq_id}->{SEQ}]));
            }
            close($out_fh);
            $cluster_count++;
            $large_clusters++;
        } 
        else 
        {
            foreach my $seq_id (@array) 
            {
                unless (exists $seqH->{$seq_id}) 
                {
                    print STDERR "Warning: Sequence ID $seq_id not found in mmseq_all_seqs.fasta\n";
                    next;
                }
                &gjoseqlib::print_alignment_as_fasta($out2_fh, ([$seq_id, $seqH->{$seq_id}->{ANNO}, $seqH->{$seq_id}->{SEQ}]));
                $leftover_seqs++;
            }
            $small_clusters++;
        }
    }
    close($out2_fh);
    
    print STDERR "Created $large_clusters cluster files\n";
    print STDERR "Put $leftover_seqs sequences from $small_clusters small clusters into Leftover_Seqs.aa\n";
    
    # Only align leftover sequences if there are any
    if ($leftover_seqs > 0 && -s "Leftover_Seqs.aa") 
    {
        print STDERR "Aligning leftover sequences with MAFFT...\n";
        my $ret = system "mafft --thread 24 --quiet --reorder Leftover_Seqs.aa > Leftover_Seqs.fa";
        
        if ($ret != 0) 
        {
            die "MAFFT execution failed with return code $ret: $!\n";
        }
        
        print STDERR "MAFFT alignment of leftover sequences completed successfully\n";
    }
    else
    {
        print STDERR "No leftover sequences to align, skipping MAFFT step\n";
        # Create an empty Leftover_Seqs.fa file for consistency
        open(my $empty_fh, ">Leftover_Seqs.fa") or die "Cannot create empty Leftover_Seqs.fa: $!\n";
        close($empty_fh);
    }
    
    # Create necessary directories
    mkdir("alis") unless -d "alis";
    mkdir("corrected_alis") unless -d "corrected_alis";
    mkdir("pssms") unless -d "pssms";
    
    print STDERR "Created output directories: alis, corrected_alis, pssms\n";
    
    return 1;
}

1; # End of module