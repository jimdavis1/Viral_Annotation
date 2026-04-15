#! /usr/bin/env perl
use strict;
use JSON::XS;
use File::Slurp;
use Data::Dumper;
use Getopt::Long;
use Cwd;
use gjoseqlib;
use Proc::ParallelLoop;
use GenomeTypeObject;


my $usage = 'annotate_viral_taxon.pl [options]

    Annotate viral genomes from BVBRC or from a local contigs directory.

    Modes:
      BVBRC mode    : requires -i Taxon_Name
      Contig dir mode: requires -contigs Contigs_Dir and -meta Metadata_File
                       -i is optional in this mode; if supplied, min/max genome
                       length will be read from the JSON file and used to filter
                       contig files before annotation.

    General options:
      -h        Help
      -t        Threads (D = 24)
      -j        Full path to the JSON opts file
                (D = /home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json)

    Output directory options:
      -c        Contig GTO folder (D = Contig_GTOs)
      -a        Annotation GTO folder (D = Anno_GTOs)
      -q        Quality GTO folder (D = Quality_GTOs)
      -d        Contigs download folder, BVBRC mode only (D = Contigs)

    BVBRC mode options:
      -i        Taxon name to query BVBRC (also usable in contig dir mode
                to pull min/max lengths from the JSON file)
      -exclude  File of genome IDs to exclude
      -min      Force a min genome length integer (must pair with -max)
      -max      Force a max genome length integer (must pair with -min)

      NOTE: BVBRC mode filters genomes by the length of a single BVBRC genome
      record.  Segmented viruses whose segments are stored as separate genome
      records in BVBRC will NOT be merged; each segment will be treated as an
      independent genome.  Use contig dir mode (-contigs/-meta) if you need to
      annotate pre-assembled multi-segment genomes as a single unit.

    Contig dir mode options:
      -contigs  Path to the directory of contig files
      -meta     Metadata file with columns: file_name TAB genome_name TAB taxon_id

    Length filtering (contig dir mode):
      If -min/-max are supplied explicitly, or if -i is supplied and the taxon
      has segments defined in the JSON, contig files whose total sequence length
      falls outside [min, max] are skipped and written to a rejected file.

    Output layout:
      Contig_GTOs/   - one GTO per genome created from raw contigs
      Anno_GTOs/     - annotated GTOs (PSSM + transcript edit + splice variants),
                       plus per-genome .stdout.txt and .stderr.txt from the PSSM step
      Quality_GTOs/  - quality-assessed GTOs plus .feature_quality and .contig_quality files

    Notes:
      The annotation pipeline pipes annotate_by_viral_pssm-GTO.pl into
      get_transcript_edited_features.pl and get_splice_variant_features.pl.
      Both downstream programs skip their processing automatically when the
      viral family has no transcript-edited or splice-variant features.

';


my $contig_gto  = "Contig_GTOs";
my $anno_gto    = "Anno_GTOs";
my $quality_gto = "Quality_GTOs";
my $contigD     = "Contigs";
my $threads     = 24;
my $json        = "/home/jjdavis/bin/Viral_Annotation/Viral_PSSM.json";
my $base        = getcwd;

my ($help, $taxon, $min, $max, $exF, $contigDir, $metaF);

my $opts = GetOptions(
	'h'         => \$help,
	'i=s'       => \$taxon,
	't=i'       => \$threads,
	'j=s'       => \$json,
	'a=s'       => \$anno_gto,
	'q=s'       => \$quality_gto,
	'c=s'       => \$contig_gto,
	'd=s'       => \$contigD,
	'contigs=s' => \$contigDir,
	'meta=s'    => \$metaF,
	'min=i'     => \$min,
	'max=i'     => \$max,
	'exclude=s' => \$exF,
);

if ($help)
{
	die "$usage\n";
}

# Validate mode
if ($contigDir || $metaF)
{
	unless ($contigDir && $metaF)
	{
		die "Contig dir mode requires both -contigs and -meta\n\n$usage\n";
	}
}
else
{
	unless ($taxon)
	{
		die "Must supply -i Taxon_Name for BVBRC mode, or use -contigs/-meta for contig dir mode\n\n$usage\n";
	}
}

if (($min && !$max) || ($max && !$min))
{
	die "Must declare both -min and -max\n\n$usage\n";
}


# Load the exclude hash if provided
my %exclude;
if ($exF)
{
	open(IN, "<$exF") or die "Could not open exclude file: $exF\n";
	while (<IN>)
	{
		chomp;
		$exclude{$_}++;
	}
	close IN;
}


# Read the JSON options file
open(IN, "<$json") or die "Cannot find JSON options file: $json\n";
my $options = decode_json(scalar read_file(\*IN));
close IN;


# Determine min/max genome length from JSON when not supplied on the command line.
# This applies to both modes when -i is provided.
unless ($min && $max)
{
	if ($taxon)
	{
		unless (exists $options->{$taxon}->{segments})
		{
			die "No segments entry for '$taxon' in JSON\n";
		}
		$max = 0;
		$min = 0;
		foreach my $segment (values %{$options->{$taxon}->{"segments"}})
		{
			$max += $segment->{"max_len"};
			$min += $segment->{"min_len"};
		}
	}
}


# Create output directories
mkdir $contig_gto;
mkdir $anno_gto;
mkdir $quality_gto;


# Collect the list of items to process and their metadata.
# In BVBRC mode  : keys are genome IDs, metadata comes from the BVBRC query.
# In contig mode : keys are filenames, metadata comes from the -meta file.

my @items;          # genome IDs (BVBRC) or contig filenames (contig dir)
my %genome_meta;    # item => { NAME => ..., TAX => ... }


if ($contigDir)
{
	# -------------------------------------------------------------------------
	# Contig directory mode: read metadata file, then optionally filter by
	# total contig length using gjoseqlib::read_fasta.
	# -------------------------------------------------------------------------

	my $label = $taxon // "contigs";
	open(GOOD, ">$label.processed") or die "Cannot open $label.processed\n";
	open(BAD,  ">$label.rejected")  or die "Cannot open $label.rejected\n";

	open(IN, "<$metaF") or die "Could not open metadata file: $metaF\n";
	while (<IN>)
	{
		chomp;
		my ($file, $name, $tax) = split /\t/;

		next if exists $exclude{$file};

		# Sum the lengths of all contigs in the file
		my @seqs     = gjoseqlib::read_fasta("$contigDir/$file");
		my $total_len = 0;
		foreach my $seq (@seqs)
		{
			$total_len += length($seq->[2]);
		}

		# Apply length filter if min/max are defined
		if ($min && $max && (($total_len < $min) || ($total_len > $max)))
		{
			print BAD "$file\t$total_len\t$name\n";
			next;
		}

		push @items, $file;
		$genome_meta{$file}->{NAME} = $name;
		$genome_meta{$file}->{TAX}  = $tax;
		print GOOD "$file\t$total_len\t$name\n";
	}
	close IN;
	close GOOD;
	close BAD;
}
else
{
	# -------------------------------------------------------------------------
	# BVBRC download mode: query BVBRC for genomes of the right length
	# -------------------------------------------------------------------------
	mkdir $contigD;

	open(IN,   "echo $taxon | query_PATRIC_bob.pl -c genome -i taxon_lineage_names -r \"genome_id genome_length genome_name\" | ");
	open(GOOD, ">$taxon.processed") or die "Cannot open $taxon.processed\n";
	open(BAD,  ">$taxon.rejected")  or die "Cannot open $taxon.rejected\n";

	while (<IN>)
	{
		chomp;
		my ($id, $len, $name) = split /\t/;
		print "$_\n";

		unless (exists $exclude{$id})
		{
			if (($len < $max) && ($len > $min))
			{
				push @items, $id;
				$genome_meta{$id}->{NAME} = $name;
				$genome_meta{$id}->{TAX}  = $id;   # taxon derived from genome_id below
				print GOOD "$id\t$len\t$name\n";
			}
			else
			{
				print BAD "$id\t$len\t$name\n";
			}
		}
	}
	close IN;
	close GOOD;
	close BAD;
}


# ============================================================================
# Main parallel processing loop
# ============================================================================

pareach(
	\@items,

	sub
	{
		my $item = shift @_;

		my $name = $genome_meta{$item}->{NAME};
		my $tax  = $genome_meta{$item}->{TAX};
		$tax =~ s/\..+//g;   # strip version suffix to get integer taxon ID

		# Derive a clean ID: strip common contig-file extensions if present
		my $id = $item;
		$id =~ s/\.(?:dna|fasta|fa|fna|contigs?)$//;


		# ------------------------------------------------------------------
		# Step 1: Build the contig GTO from raw sequences
		# ------------------------------------------------------------------
		if ($contigDir)
		{
			system "rast-create-genome "
			     . "--ncbi-taxonomy-id $tax "
			     . "--scientific-name \"$name\" "
			     . "--domain Viruses "
			     . "--genetic-code 11 "
			     . "--contigs $base/$contigDir/$item "
			     . "> $base/$contig_gto/$id.contig.gto";
		}
		else
		{
			system "echo $item | BVBRC_clean_contigs.pl -d $base/$contigD";
			system "rast-create-genome "
			     . "--ncbi-taxonomy-id $tax "
			     . "--scientific-name \"$name\" "
			     . "--domain Viruses "
			     . "--genetic-code 11 "
			     . "--contigs $base/$contigD/$item.contigs "
			     . "> $base/$contig_gto/$id.contig.gto";
		}


		# ------------------------------------------------------------------
		# Step 2: Annotate and post-process
		#
		# Pipeline:
		#   annotate_by_viral_pssm-GTO.pl  (PSSM-based feature calls)
		#     -> get_transcript_edited_features.pl  (skips if family has none)
		#     -> get_splice_variant_features.pl     (skips if family has none)
		#
		# The final GTO from this pipeline is the Anno GTO.
		# annotate_by_viral_pssm-GTO.pl writes $id.stdout.txt and
		# $id.stderr.txt to the working directory; we move those into
		# Anno_GTOs for traceability.
		# ------------------------------------------------------------------

		system "annotate_by_viral_pssm-GTO.pl -x $id "
		     . "< $base/$contig_gto/$id.contig.gto "
		     . "| get_transcript_edited_features.pl "
		     . "| get_splice_variant_features.pl "
		     . "> $base/$anno_gto/$id.anno.gto";

		system "mv $base/$id.stdout.txt $base/$anno_gto/" if -f "$base/$id.stdout.txt";
		system "mv $base/$id.stderr.txt $base/$anno_gto/" if -f "$base/$id.stderr.txt";


		# ------------------------------------------------------------------
		# Step 3: Quality assessment
		#
		# Reads the finished Anno GTO and writes a Quality GTO plus two
		# tabular report files that we move into Quality_GTOs.
		# ------------------------------------------------------------------

		system "viral_genome_quality.pl "
		     . "-p $id "
		     . "-i $base/$anno_gto/$id.anno.gto "
		     . "-o $base/$quality_gto/$id.qual.gto";

		system "mv $base/$id.feature_quality $base/$quality_gto/" if -f "$base/$id.feature_quality";
		system "mv $base/$id.contig_quality  $base/$quality_gto/" if -f "$base/$id.contig_quality";
	},

	{ Max_Workers => $threads }
);