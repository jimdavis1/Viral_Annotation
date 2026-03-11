#! /usr/bin/env perl
use strict;
use Data::Dumper;
use Getopt::Long;
use GenomeTypeObject;
use JSON::XS;
use gjoseqlib;

my $usage = 'gto_to_fasta.pl  [options] < genome.gto

    -h help
    -f filter string (case-insensitive, partial match against annotation)
    -g good quality genomes only (prints genome id and quality to stderr if Poor)
    -q good quality features only (prints feature id, annotation, and quality to stderr if Poor)
    -u exclude features with undefined quality (prints feature id, annotation, and quality to stderr)
    -d return DNA sequences instead of proteins

    Reads a single GTO from stdin and writes protein or DNA features
    as FASTA to stdout.

';

my ($help, $filter, $good_only, $good_features, $exclude_undef, $dna_mode);
my $opts = GetOptions( 'h'   => \$help,
                       'f=s' => \$filter,
                       'g'   => \$good_only,
                       'q'   => \$good_features,
                       'u'   => \$exclude_undef,
                       'd'   => \$dna_mode ) or die "$usage\n";

if ($help) { die "$usage\n"; }

my $genome = GenomeTypeObject->create_from_file(\*STDIN);
$genome or die "Error reading and parsing GTO from stdin\n";

my $gid  = $genome->{"id"};
my $qual = $genome->{"quality"}->{"genome_quality"};

if ($good_only && $qual !~ /Good/)
{
    print STDERR "$gid\t$qual\n";
    exit 0;
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

for my $i (0..$#features)
{
    my $type  = $features[$i]->{"type"};
    my $id    = $features[$i]->{"id"};
    my $anno  = $features[$i]->{"function"};
    my $fqual = $features[$i]->{"feature_quality"};

    if (($type =~ /CDS/) || ($type =~ /mat_peptide/))
    {
        next if ($filter && (!defined $anno || $anno !~ /$filter/i));

        if (!defined $fqual)
        {
            print STDERR "$id\t$anno\tFeature Quality Undefined\n";
            next if $exclude_undef;
        }
        elsif ($good_features && $fqual !~ /Good/)
        {
            print STDERR "$id\t$anno\tPoor\n";
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
                    $seq .= &revcomp($segment);
                }
            }
        }
        else
        {
            $seq = $features[$i]->{"protein_translation"};
            next unless $seq;
            $seq =~ s/\*$//g;
        }

        push @seqs, ([$id, $anno, $seq]);
    }
}

&gjoseqlib::print_alignment_as_fasta(\*STDOUT, @seqs) if @seqs;


sub revcomp
{
    my ($seq) = @_;
    $seq = reverse $seq;
    $seq =~ tr/ACGTacgt/TGCAtgca/;
    return $seq;
}