# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

This repository contains tools for viral genome annotation using Position-Specific Scoring Matrices (PSSMs). It covers Bunyavirales, Coronaviridae, Filoviridae, Orthomyxoviridae, Paramyxoviridae, and Pneumoviridae. The system identifies known viral proteins by BLAST-matching against curated PSSMs—it is not a de novo ORF discovery tool.

## Running the Annotation Pipeline

### Basic annotation from FASTA:
```bash
annotate_by_viral_pssm.pl -i contigs.fasta -p output_prefix
```

### Full pipeline with GTO (Genome Type Object) input/output:
```bash
# Step 1: Initial annotation
annotate_by_viral_pssm-GTO.pl -x prefix -i input.gto -o step1.gto

# Step 2: Transcript editing (for Paramyxoviridae phosphoproteins, Filoviridae glycoproteins)
get_transcript_edited_features.pl -i step1.gto -o step2.gto

# Step 3: Splice variants (for Influenza)
get_splice_variant_features.pl -i step2.gto -o step3.gto

# Step 4: Quality assessment
viral_genome_quality.pl -i step3.gto -o final.gto -p quality_report
```

### Key options for annotate_by_viral_pssm.pl:
- `-j` JSON options file (default: Viral_PSSM.json)
- `-c` Representative contigs directory (default: Viral-Rep-Contigs/)
- `-pssm` PSSMs directory (default: Viral-PSSMs/)
- `-min`/`-max` Contig length bounds (default: 300/35000)
- `-mcb` Minimum contig bitscore (default: 150)
- `-threads` BLAST threads
- Debug: `-tmp` (keep temp), `-no` with `-aa`/`-dna`/`-tbl` for stdout output

## Architecture

### Annotation Strategy
1. **Taxon identification**: BLASTn input contigs against `Viral-Rep-Contigs/` to find closest genus
2. **Feature calling**: tBLASTn each genus-specific PSSM from `Viral-PSSMs/{genus}.pssms/` against contigs
3. **Special processing**: Apply rules from `Viral_PSSM.json` (extensions, coverage cutoffs, copy numbers)
4. **Post-processing**: Transcript editing and splice variant detection for specific taxa

### Key Data Files
- `Viral_PSSM.json`: Feature parameters (bit cutoffs, coverage, copy numbers, segment info, PMIDs)
- `Viral-Rep-Contigs/`: Representative genomes per genus for initial classification
- `Viral-PSSMs/`: Hand-curated PSSMs organized as `{Genus}.pssms/` directories
- `Transcript-Editing/`: Post-edited transcript sequences for phosphoproteins/glycoproteins
- `Splice-Variants/`: Curated splice site sequences with header format `>ID SD:start-end;nt SA:start-end;nt`

### PSSM Generation (Other_Scripts/)
The `fasta-cluster-pssm-2.pl` script builds PSSMs from protein FASTA:
```bash
cat proteins.faa | fasta-cluster-pssm-2.pl -a "Annotation" -p prefix
```
Uses MMSeqs2 clustering (~80% identity), MAFFT alignment, then PSI-BLAST PSSM generation. Output: `clusters/`, `alis/`, `corrected_alis/`, `pssms/` directories.

## Dependencies

- Perl 5.38+ with: JSON::XS, File::Slurp, IPC::Run, Getopt::Long, Getopt::Long::Descriptive
- `gjoseqlib.pm` from https://github.com/TheSEED/seed_gjo/
- BV-BRC modules: GenomeTypeObject.pm, P3DataAPI
- BLAST+ 2.13.0 (blastn, tblastn)—JSON output format is version-specific
- MMSeqs2, MAFFT (for PSSM generation)

For BV-BRC internal users: `source /vol/patric3/cli/ubuntu-cli/user-env.sh`

## Environment Variable

`LOWVAN_DATA_DIR`: Override default data directory paths (default: `/home/jjdavis/bin/Viral_Annotation`)
