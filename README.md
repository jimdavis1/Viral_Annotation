# Viral Annotation
This repo contains code and data for improving viral annotation.  It currently covers the  *Paramyxoviridae*, *Bunyavirales*, and *Filoviridae*. The overall goal is to create a low-tech solution for calling viral proteins across entire viral families, and to cover cases where we do not have bespoke species-specific annotations from VIGOR4.<br>

This program is not intended to be used as a *de novo* protein or ORF discovery tool.  It is designed to find proteins that we already know to exist.  I will explain more about how it works below.

## Covered viral taxa
Aquaparamyxovirus<br>
Arenaviridae<br>
Ferlavirus<br>
Fimoviridae<br>
Hantaviridae<br>
Henipavirus<br>
Jeilongvirus<br>
Metaavulavirus<br>
Morbillivirus<br>
Nairoviridae<br>
Orthoavulavirus<br>
Orthoebolavirus<br>
Orthomarburgvirus<br>
Orthorubulavirus<br>
Peribunyaviridae<br>
Paraavulavirus<br>
Pararubulavirus<br>
Phasmaviridae<br>
Phenuiviridae<br>
Respirovirus<br>
Tospoviridae<br>

## Dependencies

Unless otherwise stated, the programs described in this repo are written and tested in in perl (v5.38.0).

The script(s) have the following dependencies:<br>
External CPAN perl modules:<br>

Data::Dumper<br>
File::Copy<br>
File::Path<br>
File::SearchPath<br>
File::Slurp<br>
Getopt::Long<br>
IPC::Run<br>
JSON::XS<br>
Time::HiRes<br>
<br>

It uses `gjoseqlib.pm` which is perl module that was written by Gary Olsen at the University of Illinois.  This is used for basic sequence manipulation.  You can get the latest version of this module by downloading it from Gary's repo here: https://github.com/TheSEED/seed_gjo/  <br>

There are 2 BV-BRC supported modules as well:<br>
GenomeTypeObject.pm (https://github.com/BV-BRC/p3_core/blob/master/lib/GenomeTypeObject.pm).  This perl module contains all the necessary tooling to read and write to and from a Genome Type Object, which is a JSON representation of a genome and all analysis events.<br>

The P3DataAPI "https://www.bv-brc.org/docs/cli\_tutorial/command\_list/P3DataAPI.html) is used to generate unique ids for features in the BVBRC annotation system.<br>


The program(s) run the blast suite of tools from NCBI.  The current version requires:<br>
`blastn: 2.13.0+`
`tblastn: 2.13.0+`<br>

It has not been tested on other versions of BLAST.  Note that it uses the JSON output of BLAST and other versions of BLAST have slightly different JSON formatting.  <br>

For internal users, source:<br>
`/vol/patric3/cli/ubuntu-cli/user-env.sh`<br>

## Repo Contents
`annotate_by_viral_pssm.pl` the perl script that runs the BLASTs and calls the ordinary proteins (CDSs), mature peptides, and location-based features such as RNAs. <br>
`annotate_by_viral_pssm-GTO.pl` this perl script runs annotate_by_viral_pssm.pl and creates a GTO as output. Note that its help options are slightly different.<br>
It is run in the following way: `annotate_by_viral_pssm-GTO.pl  -x [file_prefix] -i Input.gto -o Output.gto` with other options in the help menu.<br>
`get_transcript_edited_features.pl` This script reads a GTO and finds sequences that have undergone transcript editing, updating the resulting CDS or mat_peptide.<br>
`viral_genome_quality.pl`  This script reads the GTO and evaluates the genome quality based on CDSs and mat_peptide features present, and their copy number.  It also evaluates the contigs based on copy number and the proteins they encode.<br>
`Viral_PSSM.json`  This file contains BLAST and ORF parameters per feature.<br>
`Viral-Rep-Contigs` This is the directory of representative contigs that guides the program to the closest set of PSSMs.<br>
`Viral-PSSMs` the directory of hand curated PSSMS per family or genus. There may be more than one PSSM per protein.<br>
`Viral-Alignments`  This directory contains the alignments that correspond to each PSSM.  This is not used by the program, but it is useful for keeping track of the source data used to build each PSSMs.<br>
`Other-Scripts` is a directory other non-essential but useful scipts and files related to the development and management of this tool.  It currently contains a program called, `list_annos_from_pssms.pl` which will dump the annotation for each pssm.<br><br>

## How to run annotate_by_viral_pssm.pl
`annotate_by_viral_pssm.pl [options] -i subject_contig(s).fasta`<br><br>

Options include:
```

		-h help
		-i input subject contigs in fasta format
		-t declare a temp file (d = random)
		-tax declare a taxonomy id (D = 11158 )
		-g Genome name (D = Paramyxoviridae);
		-k Keep internal stop codons (D = off) if you think that your genome will have stops
		   within the PSSM, but still want to make a call over that region creating a pseudo gene.

		-min minimum contig	length (d = 1000) # otherwise the genome is rejected
		-max maximum contig length (d = 25000) # for reference Measles is 15894 and beilong is 19,212

        -opt Options file in JSON format which carries data for match (D = /home/jjdavis/Viral_PSSM.json)
		-l Representative contigs directory (D = /home/jjdavis/bin/Viral-Rep-Contigs)
		-p Base directory of PSSMs   (D = /home/jjdavis/bin/Viral-PSSMs)
	      Note that this is set up as a directory of pssms
	      right now this is hardcoded as: "virus".pssms within this directory.
```
Hard-coded locations currently exist as the defaults for -opt, -l, and -p.  Since that is annoying, you 
might want to run something like:<br>
`perl -i -pe 's/\/home\/jjdavis\/bin/the path to your bin/g' annotate_by_viral_pssm.pl`, or you could edit lines 70-72 by hand (but note that these are the line numbers at the time I wrote this).<br><br>




There is also a set of debugging parameters that I use frequently:
```
 -tmp keep temp dir
 -no no output files generated, to be used in conjunction with one of the following:
  -dna print only genes to STDOUT 
  -aa print proteins to STDOUT
  -tbl print only feature table to STDOUT
  -ctbl [file name] concatenate table results to a file (for use with many genomes)
```


## How it works

The code is currently designed to work on the *Paramyxoviridae*, *Bunyavirales*, and Filoviridae, although more families are planned.  As depicted in the image below, it first performs a BLASTn against a small set of representative genomes for each genus.  Then it sorts the results by bit score and chooses the best match.<br><br>
For each genus, there is a directory of PSSMs corresponding to each known protein for that genus. The PSSMs are derived from a set of hand curated alignments. In the next step, it cycles through each directory of PSSMs (there may be more than one PSSM per protein), choosing the best tBLASTn match per pssm. <br>

Note that it assumes your genome will have the same set of proteins as the nearest match. This is why it is not intended to be used as a discovery tool.  In the event that a new protein is found, a new PSSM must be added to the PSSM directory.  <br><br>

![Anno-Strategy](https://github.com/jimdavis1/Viral_Annotation/assets/7661533/0d6a3a44-47af-40bf-852d-5ddda250ad94)

<br><br>Finally it performs any special rules on the proteins/ORFs.  These rules are currently encoded in a JSON file called `Viral_PSSM.json`. The following is a description fo the current JSON strucutre. <br>

The perl code reads the JSON file into a hashref, which has this general structure: <br><br>

```
  "Arenaviridae": {
    "segments": {
      "Small RNA Segment": {
        "max_len": 3741,
        "min_len": 3061,
        "replicon_geometry": "linear"
      },
      "Large RNA Segment": {
        "max_len": 8014,
        "min_len": 6556,
        "replicon_geometry": "linear"
      }
    "features": {
      "GPC": {
        "anno": "Pre-glycoprotein polyprotein GP complex (GPC protein)",
        "bit_cutoff": 100,
        "copy_num": 1,
        "coverage_cutoff": 0.65,
        "downstream_ext": 1,
        "feature_type": "CDS",
        "max_len": 617,
        "min_len": 455,
        "non_pssm_partner": ["Small Segment Stemloop"],
        "segment": "Small RNA Segment",
        "upstream_ext": 1
      },
      ...
```

The Viral_PSSM.json file is in a regular state of development, so this may change slightly, but the above shows an example for, *Arenaviridae*, and a single protein, GPC.  The two highest level keys are `segments`, which contains information on segments that are used for genome quality evaluation.  `features` currently contains information on CDS, mat_peptide, and RNA features. <br>

The following is a non-exhaustive description of fields that are used in the JSON<br>
`max_len and min_len` maximum or minimum length of a contig or feature that is expected and evaluated by the genome quality checker (not all features or contigs will have this)<br>
`replicon_geometry` this is not currently used, but carries info on the genometry of the replicon and inserted into the GTO by the quality tool<br>
`copy_num` expected copy number of a feature<br>
`coverage_cutoff` blast subject coverage for calling a feature<br>
`upstream_ext and downstream_ext` tells the program if it can look upstream for a met start or downstream for a stop codon<br>
`feature_type` currently CDS, mat_peptide, or RNA <br>
`segment` which segment a feature belongs to (used by quality tool)<br>
`non_pssm_partner` used for placing a location based feature.<br><br>

There are other fields that are not depicted in the example, including:<br>
`PMID` which contains the pubmed ID for one or more DLITS.  These are examples of important papers that either define the function or sequence of a feature. <br>
`"special": "transcript_edit"`  This field tells the program that an external program is being used to make a call.  In this case, `transcript_edit` is used to denote a feature that undergoes transcript editing and is found by using `get_transcript_edited_features.pl`.<br>


## General notes on the curation and development of PSSMs
### Paramyxoviridae

I have recently updated the way transcript-edited features are called by adding `get_transcript_edited_features.pl`.  This is up-to-date and evaluated for the glycoproteins of ebola, but I need to go back and update the phosphoproteins in the paramyxos.  They are currently called by splicing two ORFS, which turned out to be problematic in a few cases. 


<br>


### Tospoviridae:
I was unable to find any acceptable publications that unambiguously define the coordinates of Gn and Gc.<br>

### Fimoviridae:
I also could not find any publications clearly showing Gc and Gn.<br>

The Fimoviridae are the most poorly characterized family that I have encountered so far.  They are  multi-segmented and variable in their smaller segments. Proteins from these segments including P5, 6, 6a, 6b, 7, and 27 are all essentially uncharacterized.  They are numbered based on appearance in the genome in which they are described, but their ordering may or may not hold up as more are sequenced.  Furthermore, the proteins that have been called P5 and P6 have little to no similarity amongst themselves (usually < 35% identity) and could all have different functions in their own right.  I chose to split these into individual sets of pssms with the annotation "Fimoviridae uncharacterized protein."  We can hang an annotation on each when we learn what it does.  It is worth noting that due to the infrequency of these proteins, there are many low-occurrence uncharacterized proteins that did not get PSSMs and are not getting called.   The "P5" protein of Raspberry leaf blotch emaravirus is a good example here (fig|1980431.35.CDS.1).<br> 

In this family the quality checker will look for Segments 1-4 only, which correspond to the individual proteins L, GPC, N, and MOV, respectively.  Their segment lengths are highly variable, so the lenght cutoffs for segments 1-4 are based on the the lower length limit of the corresponding protein, and (the longest allowable gene + 0.5 X longest allowable gene) (this is arbitrary and could  be tuned).<br>

## Phasmaviridae
These are mostly insect virueses.  The set of genomes is highly diverse,with few representatives in each genus, so the pssms only represent a fraction of the true diversity.  There were a considerable number of proteins that I could not get to cluster at 50% identity. I am currently dissatisfied with this family, so as more exemplars come in, this set should eventually get recomputed. 








