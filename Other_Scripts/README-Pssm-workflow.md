# PSSM Generation Pipeline

This document explains the workflow of the `fasta-cluster-pssm-2.pl` program.  The program is an automated approach to generating Position-Specific Scoring Matrices (PSSMs) for the LowVan annotation system. It reads in FASTA formatted files of proteins that are all expected to encode isofunctional homologs and returns curated alignments and PSSMs. 

The objective of the pipeline is to create a set of PSSMs with well conserved N- and C-terminal regions, and that will result in a single BLAST HSP when searched against a new genome.  To achieve this, the pipeline conducts a series of alignment editing steps that make the alignments as compact as possible with out destroying too much of the natural variation in the sequences.  You may want to take this into consideration if you intend to use this script as an alignment or PSSM generator for other applications.  The proteins are also clustered at ~80% sequence identity, so they are probably not ideal for long distance sequence similarity searches.


## Table of Contents

1. [Overview](#overview)
2. [Pipeline Workflow](#pipeline-workflow)
3. [Input and Output](#input-and-output)
4. [Module Descriptions](#module-descriptions)
5. [Detailed Process Steps](#detailed-process-steps)
6. [N-terminal Conservation Analysis](#n-terminal-conservation-analysis)
7. [PSSM Generation](#pssm-generation)
8. [Command Line Options](#command-line-options)
9. [Directory Structure](#directory-structure)

## Overview

The FASTA Cluster PSSM pipeline is designed to process protein sequences, cluster them based on sequence similarity, generate multiple sequence alignments for each cluster, and create PSSMs.

## Pipeline Workflow

The pipeline follows these main steps:

1. Parse command line arguments
2. Read protein sequences from standard input
3. Process and filter sequences
4. Run MMSeqs clustering
5. Align each cluster with MAFFT
6. Process alignments (remove low occupancy columns, trim ends)
7. Create PSSMs for each processed alignment
8. Evaluate N-terminal conservation and recluster if necessary
9. Generate final PSSMs for reclustered alignments

## Input and Output

### Input
- Protein sequences in FASTA format (from standard input)
- Command line parameters

### Output
- Cluster FASTA files (`clusters/` directory)
- Original alignments (`alis/` directory)
- Processed alignments (`corrected_alis/` directory)
- Position-Specific Scoring Matrices (`pssms/` directory)
- Reclustered alignments (`reclustered_alis/` directory if N-terminal reclustering is performed)
- Log files (`Curation_Report`, `nterm_evaluation.log`, `nterm_reclustering.log`)

## Module Descriptions

The pipeline uses four perl modules that have to be installed in the path:

1. **FCP_Main_Utils**: Handles command line parsing, sequence processing, and MMSeqs clustering
2. **FCP_Ali_Utils**: Processes alignments, handles column filtering and PSSM creation
3. **FCP_PSSM_Utils**: Creates and processes PSSMs using BlastInterface
4. **FCP_Nterm_Utils**: Evaluates N-terminal conservation and performs reclustering

## Detailed Process Steps

### 1. Sequence Processing (`process_sequences` in FCP_Main_Utils)

This step processes the input sequences based on user options:

- Identical sequences are removed unless `-i` option is used
- Sequences with ambiguous residues are removed unless `-x` option is used
- Sequences are sorted by length (longest first)
- Sequences are written to a temporary FASTA file

### 2. MMSeqs Clustering (`run_mmseqs` in FCP_Main_Utils)

The MMSeqs2 tool is used to cluster sequences based on sequence similarity:

- Runs `mmseqs easy-cluster` with user-specified identity threshold (`-mi`) and coverage (`-mc`)
- Clusters are processed based on size
- Clusters with more sequences than the minimum (`-m`)  are retained for alignment
- Smaller clusters' sequences are collected into "Leftover_Seqs.aa"

### 3. Alignment Processing (`process_cluster_alignments` in FCP_Ali_Utils)

For each cluster:

- Sequences are aligned using MAFFT
- Original alignments are saved to the `alis/` directory
- Alignments are processed to remove:
  - Sequences with too many gaps (controlled by `-fd` option)
  - Columns with low occupancy (controlled by `-f` option)
  - Poorly conserved N-terminal columns (controlled by `-n` option)
  - Poorly conserved C-terminal columns (controlled by `-c` option)
  - Sequences starting with gaps (unless `-d` option is used)
- Processed alignments are saved to the `corrected_alis/` directory
- PSSMs are created for alignments with enough sequences

### 4. N-terminal Conservation Evaluation (`evaluate_nterm_conservation` in FCP_Nterm_Utils)

If enabled (default), this step:

- Analyzes the first 25 amino acids of each sequence in each alignment
- Identifies columns with poor conservation (below threshold set by `-p_nterm`)
- If enough poorly conserved columns are found (controlled by `-n_nterm`), the alignment is marked for reclustering

### 5. Reclustering (`recluster_by_nterm` in FCP_Nterm_Utils)

For alignments with poor N-terminal conservation:

- N-terminal regions are extracted (first 25aas)
- MMSeqs is run on these regions with a different identity threshold (`-nterm-mmseq-id`)
- New clusters are created based on N-terminal similarity
- New alignments are generated with MAFFT
- New PSSMs are created for each new alignment
- Leftover sequences (<5 sequences per cluster) are added to the leftover seqs


## PSSM Generation

Position-Specific Scoring Matrices (PSSMs) are generated using the BlastInterface module:

1. For each processed alignment, `alignment_to_pssm` is called
2. The PSSM is created using the PSI-BLAST algorithm
3. PSSM headers are processed to improve format
4. Final PSSMs are stored in the `pssms/` directory

These PSSMs can be used with PSI-BLAST or RPS-BLAST for sensitive protein similarity searches.

## Command Line Options

The pipeline accepts numerous command line options:

```
-h                 Help
-a "string"        Annotation string (becomes the title in the PSSM)
-m number          Minimum number of sequences for building PSSM (default = 5)
-f fraction        Delete columns with less than fraction bases (default = 0.33)
-n fraction        Delete columns from N-terminus until >= fraction identity (default = 0.75)
-c fraction        Delete columns from C-terminus until >= fraction identity (default = 0)
-x                 Include non-standard amino acids
-i                 Include identical sequences
-d                 Keep sequences beginning with a dash
-e number          Dash cutoff at sequence end (default = 3)
-efo fraction      End fraction occupancy (default = 0.15)
-fd fraction       Maximum fraction of dashes allowable in entire protein (default = 0.20)
-mi fraction       MMSeqs identity (default = 0.8)
-mc fraction       MMSeqs coverage (default = 0.8)
-p "prefix"        PSSM prefix
-p_nterm percent   N-terminal conservation threshold (default = 65%)
-n_nterm number    Required poor columns for reclustering (default = 3)
-nterm-mmseq-id id N-terminal reclustering identity (default = 0.7)
-no_nterm_eval     Disable N-terminal evaluation
```

## Directory Structure

The pipeline creates and uses the following directory structure:

```
Working Directory/
├── clusters/                # Cluster FASTA files
├── alis/                    # Original alignments
├── corrected_alis/          # Processed alignments
├── pssms/                   # Position-Specific Scoring Matrices
├── alis-for-reclustering/   # Alignments that need reclustering
├── reclustering/            # Temporary files for reclustering
├── reclustered_alis/        # Reclustered alignments
├── Curation_Report          # Log of alignment processing
├── nterm_evaluation.log     # Log of N-terminal evaluation
├── nterm_reclustering.log   # Log of reclustering
├── Leftover_Seqs.aa         # Sequences from small clusters (unaligned)
└── Leftover_Seqs.fa         # Sequences from small clusters (aligned)
```

