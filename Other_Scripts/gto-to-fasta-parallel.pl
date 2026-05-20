#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
use JSON::XS;
use gjoseqlib;
use Parallel::ForkManager;
use File::Temp qw(tempfile);


my $usage = 'gto_to_fasta.pl  [options] < genome.gto
             gto_to_fasta.pl  [options] -D /path/to/gto/dir > combined.fasta

    -h      help
    -f STR  filter string (case-insensitive, partial match against annotation)
    -g      good quality genomes only (prints genome id and quality to stderr if Poor)
    -q      good quality features only (prints feature id, annotation, and quality to stderr if Poor)
    -u      exclude features with undefined quality (prints feature id, annotation, and quality to stderr)
    -d      return DNA sequences instead of proteins
    -D DIR  process all *.gto files in DIR in parallel and write combined FASTA to stdout
    -p INT  number of parallel workers when using -D (default: number of CPU cores)

    Single-genome mode: reads one GTO from stdin, writes protein or DNA features as FASTA to stdout.
    Directory mode (-D): reads every *.gto file in DIR in parallel; all sequences are merged into
    a single FASTA on stdout.  Genome/feature quality filters and all other options still apply.

';

my ($help, $filter, $good_only, $good_features, $exclude_undef, $dna_mode, $dir_mode, $parallel);

my $opts = GetOptions(
    'h'   => \$help,
    'f=s' => \$filter,
    'g'   => \$good_only,
    'q'   => \$good_features,
    'u'   => \$exclude_undef,
    'd'   => \$dna_mode,
    'D=s' => \$dir_mode,
    'p=i' => \$parallel,
) or die "$usage\n";

if ($help) { die "$usage\n"; }

# ── directory (parallel) mode ─────────────────────────────────────────────────
if ($dir_mode)
{
    -d $dir_mode or die "Directory not found: $dir_mode\n";

    my @gto_files = sort glob("$dir_mode/*.gto");
    @gto_files or die "No *.gto files found in $dir_mode\n";

    # default workers = CPU cores, fall back to 4
    unless (defined $parallel)
    {
        $parallel = cpu_count();
    }

    # Each child writes its sequences to a numbered temp file so we can
    # concatenate them in stable (sorted-filename) order afterwards.
    my $tmpdir = File::Temp->newdir(CLEANUP => 1);
    my @tmp_paths;          # parallel index -> temp file path

    for my $i (0 .. $#gto_files) {
        push @tmp_paths, sprintf("%s/%05d.fasta", $tmpdir, $i);
    }

    my $pm = Parallel::ForkManager->new($parallel);

    my $total     = scalar @gto_files;
    my $completed = 0;

    $pm->run_on_finish(sub {
        $completed++;
        printf STDERR "\rProgress: %d / %d  (%.1f%%)",
            $completed, $total, 100 * $completed / $total;
        print STDERR "\n" if $completed == $total;
    });

    for my $i (0 .. $#gto_files)
    {
        my $gto_file = $gto_files[$i];
        my $tmp_path = $tmp_paths[$i];

        $pm->start and next;   # fork

        # ── child ──────────────────────────────────────────────────────────
        open(my $fh_in,  '<', $gto_file)  or die "Cannot open $gto_file: $!\n";
        open(my $fh_out, '>', $tmp_path)  or die "Cannot open $tmp_path: $!\n";

        my $genome = GenomeTypeObject->create_from_file($fh_in);
        close $fh_in;

        if ($genome)
        {
            my @seqs = extract_seqs($genome, $gto_file, 0);
            gjoseqlib::print_alignment_as_fasta($fh_out, @seqs) if @seqs;
        }
        else
        {
            print STDERR "Warning: could not parse $gto_file — skipping\n";
        }

        close $fh_out;
        $pm->finish;
        # ── end child ──────────────────────────────────────────────────────
    }

    $pm->wait_all_children;

    # concatenate all temp files to stdout in sorted order
    for my $tmp (@tmp_paths)
    {
        if (-s $tmp)   # skip empty files (filtered-out genomes)
        {
            open(my $fh, '<', $tmp) or die "Cannot read $tmp: $!\n";
            print while <$fh>;
            close $fh;
        }
    }

    exit 0;
}

# ── single-genome stdin mode ──────────────────────────────────────────────────
my $genome = GenomeTypeObject->create_from_file(\*STDIN);
$genome or die "Error reading and parsing GTO from stdin\n";

my @seqs = extract_seqs($genome, "stdin");
gjoseqlib::print_alignment_as_fasta(\*STDOUT, @seqs) if @seqs;

exit 0;


# ─────────────────────────────────────────────────────────────────────────────
# extract_seqs($genome, $source_label)
#
# Applies all quality/filter/mode options and returns an array of
# [id, annotation, sequence] triples ready for gjoseqlib output.
# ─────────────────────────────────────────────────────────────────────────────
sub extract_seqs
{
    my ($genome, $source, $verbose) = @_;
    $verbose //= 1;

    my $gid  = $genome->{"id"};
    my $qual = $genome->{"quality"}->{"genome_quality"};

    if ($good_only && $qual !~ /Good/)
    {
        print STDERR "$gid\t$qual\n" if $verbose;
        return ();
    }

    # build contig lookup hash if we need DNA
    my %contigs;
    if ($dna_mode)
    {
        foreach my $contig (@{ $genome->{"contigs"} })
        {
            $contigs{ $contig->{"id"} } = $contig->{"dna"};
        }
    }

    my @features = @{ $genome->{"features"} };
    my @seqs;

    for my $i (0 .. $#features)
    {
        my $type  = $features[$i]->{"type"};
        my $id    = $features[$i]->{"id"};
        my $anno  = $features[$i]->{"function"};
        my $fqual = $features[$i]->{"feature_quality"};

        next unless ($type =~ /CDS/ || $type =~ /mat_peptide/);

        next if ($filter && (!defined $anno || $anno !~ /$filter/i));

        if (!defined $fqual)
        {
            print STDERR "$id\t$anno\tFeature Quality Undefined\n" if $verbose;
            next if $exclude_undef;
        }
        elsif ($good_features && $fqual !~ /Good/)
        {
            print STDERR "$id\t$anno\tPoor\n" if $verbose;
            next;
        }

        my $seq;
        if ($dna_mode)
        {
            my $locations = $features[$i]->{"location"};
            $seq = "";
            for my $loc (@$locations)
            {
                my ($contig_id, $start, $strand, $length) = @$loc;
                my $contig_seq = $contigs{$contig_id};
                unless (defined $contig_seq)
                {
                    warn "Contig $contig_id not found for feature $id\n";
                    next;
                }

                if ($strand eq "+")
                {
                    $seq .= substr($contig_seq, $start - 1, $length);
                }
                else
                {
                    my $segment = substr($contig_seq, $start - $length, $length);
                    $seq .= revcomp($segment);
                }
            }
        }
        else
        {
            $seq = $features[$i]->{"protein_translation"};
            next unless $seq;
            $seq =~ s/\*$//g;
        }

        push @seqs, [$id, $anno, $seq];
    }

    return @seqs;
}


sub cpu_count
{
    # try nproc (Linux coreutils)
    my $n = do { local $_ = `nproc 2>/dev/null`; chomp; $_ + 0 };
    return $n if $n > 0;

    # try counting /proc/cpuinfo processor entries
    if (open my $fh, '<', '/proc/cpuinfo') {
        $n = scalar grep { /^processor\s*:/ } <$fh>;
        close $fh;
        return $n if $n > 0;
    }

    return 8;   # sensible default
}

sub revcomp
{
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}