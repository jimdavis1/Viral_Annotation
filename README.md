# Viral Annotation
This repo contains code and data for improving viral annotation.  It currently covers the  *Paramyxoviridae*, *Bunyavirales*, and *Filoviridae*. The overall goal is to create a low-tech solution for calling viral proteins across entire viral families and to cover cases where we do not have bespoke species-specific annotations from VIGOR4.<br>

This program is not intended to be used as a *de novo* protein or ORF discovery tool.  It is designed to find proteins that we already know to exist.  


## Covered viral taxa

### Bunyavirales:
Arenaviridae<br>
Fimoviridae<br>
Hantaviridae<br>
Nairoviridae<br>
Peribunyaviridae<br>
Phasmaviridae<br>
Phenuiviridae<br>
Tospoviridae<br>

### Filoviridae:
Orthoebolavirus<br>
Orthomarburgvirus<br>

### Orthomyxoviridae:
Influenza A virus<br>

### Paramyxoviridae:
Aquaparamyxovirus<br>
Ferlavirus<br>
Henipavirus<br>
Jeilongvirus<br>
Metaavulavirus<br>
Morbillivirus<br>
Narmovirus<br>
Orthoavulavirus<br>
Orthorubulavirus<br>
Paraavulavirus<br>
Pararubulavirus<br>
Respirovirus<br>

### Pneumoviridae:
Orthopneumovirus<br>
Metapneumovirus<br>
Orthopneumovirus muris<br>


## Dependencies

Unless otherwise stated, the programs described in this repo are written and tested in in perl (v5.38.0).

The script(s) have the following dependencies:<br>

External CPAN perl modules:<br>

Data::Dumper<br>
File::(including Copy, Path, SearchPath, Temp, Slurp)<br>
Getopt(Long and Descriptive)<br>
IPC::Run<br>
JSON::XS<br>
Time::HiRes<br>
Cwd<br>

It also uses `gjoseqlib.pm` which is perl module that was written by Gary Olsen at the University of Illinois.  This is used for basic sequence manipulation.  You can get the latest version of this module by downloading it from Gary's repo here: https://github.com/TheSEED/seed_gjo/  <br>

There are 2 BV-BRC supported modules as well:<br>
GenomeTypeObject.pm (https://github.com/BV-BRC/p3_core/blob/master/lib/GenomeTypeObject.pm).  This perl module contains all the necessary tooling to read and write to and from a Genome Type Object, which is a JSON representation of a genome and all analysis events.<br>

The P3DataAPI "https://www.bv-brc.org/docs/cli\_tutorial/command\_list/P3DataAPI.html) is used to generate unique ids for features in the BV-BRC annotation system.<br>


The program(s) run the blast suite of tools from NCBI.  The current version requires:<br>
`blastn: 2.13.0+`<br>
`tblastn: 2.13.0+`<br>

It is not guaranteed to work on other versions of BLAST.  It uses the JSON output of BLAST and other versions have slightly different JSON structures.  <br>

For internal users, source:<br>
`/vol/patric3/cli/ubuntu-cli/user-env.sh`<br>

## Repo Contents
`annotate_by_viral_pssm.pl` the perl script that runs the BLASTs and calls the ordinary proteins (CDSs), mature peptides, and location-based features such as RNAs. <br><br>

`annotate_by_viral_pssm-GTO.pl` this perl script runs annotate_by_viral_pssm.pl and creates a GTO as output. Note that its help options are slightly different.<br>
It is run in the following way: `annotate_by_viral_pssm-GTO.pl  -x [file_prefix] -i Input.gto -o Output.gto` with other options in the help menu.<br><br>

`get_transcript_edited_features.pl` This script reads a GTO (that has already been processed by annotate_by_viral_pssm-GTO.pl) and finds sequences that have undergone transcript editing, updating the resulting CDS or mat_peptide to ensure we get the correct protein sequence.<br><br>

`get_splice_variant_features.pl` This script reads a GTO (that has already been processed by annotate_by_viral_pssm-GTO.pl) and finds sequences that are the result of splicing, updating the resulting CDS to ensure we get the correct protein sequence.<br><br>

`viral_genome_quality.pl`  This script reads the GTO and evaluates the genome quality based on CDSs and mat_peptide features present, and their copy number.  It also evaluates the contigs based on copy number and the proteins they encode.  It is intended to be run downstream of annotate_by_viral_pssm-GTO.pl and get_transcript_edited_features.pl.<br><br>

`Viral_PSSM.json`  This file contains BLAST and ORF calling parameters per feature.<br><br>

`Viral-Rep-Contigs` This is the directory of representative contigs that guides the program to the closest set of PSSMs.<br><br>

`Viral-PSSMs` This is the directory of hand curated PSSMS per taxon. <br><br>

`Transcript-Editing` This directory contains fasta files of hand-curated transcripts (post editing).<br><br>

`Splice-Variants` This directory contains fasta hand-curated fasta files with splice site locations.<br><br>

`PSSM-Alignments`  This directory is not used by any program, but it contains the alignments that correspond to each PSSM.<br><br>

`Other-Scripts` is a directory other non-essential but useful scipts and files related to the development and management of these tools.  It currently contains several readme files s and the tools for building pssms.<br><br>

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


## Step 1.  Calling features based on PSSMs

The code is currently designed to work on the *Paramyxoviridae*, *Bunyavirales*, *Filoviridae*, and *Pneumoviridae*, although more taxa are planned.  As depicted in the image below, it first performs a BLASTn against a small set of representative genomes for each genus.  Then it sorts the results by bit score and chooses the best match.<br>

For each genus, there is a directory of PSSMs corresponding to each known protein for that genus. The PSSMs are derived from a set of hand curated alignments. In the next step, it cycles through each directory of PSSMs (there may be more than one PSSM per protein), choosing the best tBLASTn match per pssm. <br>

Note that it assumes your genome will have the same set of proteins as the nearest match. This is why it is not intended to be used as a discovery tool.  In the event that a new protein is found, a new PSSM must be added to the PSSM directory.  <br><br>

![Anno-Strategy](https://github.com/jimdavis1/Viral_Annotation/assets/7661533/0d6a3a44-47af-40bf-852d-5ddda250ad94)

<br><br>Finally it performs any special rules on the proteins/ORFs.  These rules are currently encoded in a JSON file called `Viral_PSSM.json`. The following is a description of the current JSON strucutre.<br><br>

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

The Viral_PSSM.json file is in a regular state of development, so this may change slightly, but the above shows an example for, *Arenaviridae*, and a single protein, GPC.  The two highest level keys are `segments`, which contains information on segments that are used for genome quality evaluation and `features`, which currently contains information on CDS, mat_peptide, and RNA features. <br>

The following is a non-exhaustive description of fields that are used in the JSON<br><br>

`max_len and min_len` maximum or minimum length of a contig or feature that is expected and evaluated by the genome quality checker (not all features or contigs will have this).  The boundaries are currently very crude, and not bound by any sort of statistics, but effective<br><br>

`replicon_geometry` this is not currently used, but carries info on the genometry of the replicon and inserted into the GTO by the quality tool<br><br>

`copy_num` expected copy number of a feature<br><br>

`coverage_cutoff` blast subject coverage for calling a feature<br><br>

`upstream_ext and downstream_ext` tells the program if it can look upstream for a Met start or downstream for a stop codon<br><br>

`feature_type` currently CDS, mat_peptide, or RNA <br><br>

`segment` which segment a feature belongs to (used by the quality tool)<br><br>

`non_pssm_partner` used for placing a location based feature<br><br>

There are other fields that are not depicted in the example, including:<br>
`PMID` which contains the PubMed ID for one or more DLITS.  A DLIT is an example of an important paper that either defines the function or sequence of a feature. <br><br>

`"special": "transcript_edit"`  This field tells the program that an external program is being used to make a call.  In this case, `transcript_edit` is used to denote a feature that undergoes transcript editing and is found by using `get_transcript_edited_features.pl`.  `splice` is also a valid field and triggers the search for splice variant features.<br><br>


## Get Transcript Edited Features
Transcript editing is a phenomenon that occurs in the phosphoproteins of the *Paramyxoviridae* and the glycoproteins of the *Filoviridae*.  It occurs when the RNA-Dependent RNA polymerase encounters a region of low complexity and pauses.  The pause allows for the insertion of one or more new nucleotides into the transcript, which causes a frame shift. Thus, the amino acid sequence is not a direct translation of what is encoded in the genome.  We solve this problem by hand-curating a set of transcripts in their post-editing state. These are found in the `Transcript-Editing` directory.  We then BLAST these against the the contig, and for BLASTn matches with high enough scores, the alignment gap is filled in using the nucleotide sequence of closest curated transcript.  Currently in order to do this, the following strict BLASTn criteria must be met: <br>

1.  The match must have >= 95% nucleotide identity
2.  The match must have >= 95% query coverage
3.  The match must have <= 2 gap characters in the subject
4.  The gap characters must occur consecutively in a run <br><br>

These parameters are controlled using `--id`, `--cov`, and `--gaps` options, respectively.  The requirement for consecutive gap characters is hard-coded.
<br><br>

Because we may encounter a decent BLASTn match, but not have sufficient %identity, %query coverage, or there may be additional naturally-occurring gaps in the subject, this program will call a feature covering the BLASTn match when the above inclusion criteria are not met.  However, it will not attempt to correct the subject sequence.  Instead it will call a `partial_cds` feature and will not attempt a translation.  Parameters setting the minimum BLAST requirements for this type of feature call are `--eval`, `--lower_pid`, and `--lower_pcov`, which set the maximum BLAST e-value, and the minimum percent identity, and query coverage for consideration. <br><br>

Full usage for this program is as follows:


```	--input STR (or -i)    Input GTO
	--output STR (or -o)   Output GTO
	--cov INT (or -c)      Minimum BLASTn percent query coverage (D = 95)
	--id INT (or -p)       Minimum BLASTn percent identity  (D = 95)
	--gaps INT (or -g)     Maximum number of allowable gaps (D = 2)
	--e_val NUM (or -e)    Maximum BLASTn evalue for considering any HSP
	                       (D = 0.5)
	--lower_pid            Lower percent identity threshold for a feature
	                       call without transcript editing correction (D
	                       = 80)
	                       aka --lpi
	--lower_pcov           Lower percent query coverage for for a feature
	                       call without transcritp editing correction (D
	                       = 80)
	                       aka --lpi
	--threads INT (or -a)  Threads for the BLASTN (D = 24))
	--json STR (or -j)     Full path to the JSON opts file
	--dir STR (or -d)      Full path to the directory hand curated
	                       transcripts
	--tmp STR (or -t)      Declare name for temp dir (D = randomly named
	                       in cwd)
	--help (or -h)         Show this help message
	--debug (or -b)        Enable debugging
	
	```



## Get Splice Variant Features
The use of splicing is fairly common in viruses, and is currently necessary for calling many proteins in *Influenza*.

In order to support splice variant calling we maintain a directory of hand-curated DNA sequences with the coordinates of the splice site carefully delineated.  These are found in the `Splice-Variants` directory.  We find these these by performing a BLASTn search against the contig using our curated sequences as the query.  Then for the BLASTn matches with high enough scores, the splice is made using the curated coordinates in the fasta header.<br> 

It is worth noting that:
1.  Query sequences don't need to be be aligned, but it is easier to deal with them if they are.
2.  If query sequences are derived from an alignment they should not contain gap characters upstream of the splice junction. 
3.  All query sequences must be in the forward direction  
<br>

Fasta headers are formatted in the following way:<br>

`>valid_sequence_ID SD:SD_Region_Start-SD_Region_End;Last_nt_of_SD  SA:SA_Region_Start-SA_Region_End;First_nt_of_SA`
<br>
where SD is sequence donor, and SA is sequence acceptor.  Here is what one looks like:
`>1413195.5 SD:371-381;373 SA:491-504;503`
 <br>

Full usage for this program is as follows:

```
	--input STR (or -i)    Input GTO
	--output STR (or -o)   Output GTO
	--cov INT (or -c)      Overall Minimum BLASTn percent query coverage
	                       (D = 95)
	--id INT (or -p)       Overall Minimum BLASTn percent identity  (D =
	                       95)
	--threads INT (or -a)  Threads for the BLASTN (D = 24))
	--json STR (or -j)     Full path to the JSON opts file
	--dir STR (or -d)      Full path to the directory hand curated
	                       transcripts
	--tmp STR (or -t)      Declare name for temp dir (D = randomly named
	                       in cwd)
	--help (or -h)         Show this help message
	--debug (or -b)        Enable debugging
```

## Genome Quality Tool
As described above, the JSON file that contains information about the features also contains information about copy number of features and contigs.  The quality tool assess the the following things:<br>

1.  The number of ambigous bases per contig
2.  The number of expected segments 
3.  The length of each segment relative to what is expected
3.  The number of expected occurrences of each non-variable feature
4.  The legnth of each non-variable feature<br><br>

Currently, the tool mostly looks for CDS features, but it looks for some mat_peptides in the Filoviridae. 

The output is two tables: one is a contig report and the other is a feature report.  If any given contig or feature causes the genome quality to be "poor" the reason for the call is provided.  

Usage statement for the tool:
```viral_genome_quality.pl [-ahijop] [long options...]
	--ambiguous NUM (or -a)  Fraction of ambiguous bases, (Default = 0.01)
	--input STR (or -i)      Input GTO
	--output STR (or -o)     Output GTO
	--prefix STR (or -p)     Genome Quality File Prefix
	--json STR (or -j)       Full path to the JSON opts file
	--help (or -h)           Show this help message
```

## General remarks on the curation and development of PSSMs and the current state of the annotations
### Paramyxoviridae

I have recently updated the way transcript-edited features are called by adding `get_transcript_edited_features.pl`.  This is up-to-date and evaluated for the glycoproteins of Ebola, and the phosphoproteins in the Paramyxos.  They were originally called by splicing two BLAST HSPs, which turned out to be problematic in a few cases. DLITs that either describe the editing site, or the subsequent amino acid sequence for the transcript-edited proteins have been added to the json.  There are a handful, like Narmovirus, where I do not think protein work has been done to prove V and W, but the predicted editing site is supported by literature. At this point, all editing sites are backstopped by literature references. <br>

### Respirovirus
Note that in the respiroviruses, there is a nomenclature discrepancy regarding whether the third phosphoprotein (+2 G) is called W or D.  Currently these are all called W by the system and I have not enountered a compelling reason (other than the historical naming) to maintain the distinction between W and D.

### Tospoviridae:
I was unable to find any acceptable publications that unambiguously define the coordinates of Gn and Gc.<br>

### Fimoviridae:
I also could not find any publications clearly showing Gn and Gc.<br>

The Fimoviridae are the most poorly characterized family that I have encountered so far.  They are  multi-segmented and variable in their smaller segments. Proteins from these segments including P5, 6, 6a, 6b, 7, and 27 are all essentially uncharacterized.  They are numbered based on appearance in the genome in which they are described, but their ordering may or may not hold up as more are sequenced.  Furthermore, the proteins that have been called P5 and P6 have little to no similarity amongst themselves (usually < 35% identity) and could all have different functions in their own right.  I chose to split these into individual sets of pssms with the annotation "Fimoviridae uncharacterized protein."  We can hang an annotation on each when we learn what it does.  It is worth noting that due to the infrequency of these proteins, there are many low-occurrence uncharacterized proteins that did not get PSSMs and are not getting called.   The "P5" protein of Raspberry leaf blotch emaravirus is a good example here (fig|1980431.35.CDS.1).<br> 

In this family the quality checker will look for Segments 1-4 only, which correspond to the individual proteins L, GPC, N, and MOV, respectively.  Their segment lengths are highly variable, so the lenght cutoffs for segments 1-4 are based on the the lower length limit of the corresponding protein, and (the longest allowable gene + 0.5 X longest allowable gene) (this is arbitrary and could  be tuned).<br>

## Phasmaviridae
These are mostly insect virueses.  The set of genomes is highly diverse with few representatives in each genus, so the pssms only represent a fraction of the true diversity.  There were a considerable number of proteins that I could not get to cluster at 50% identity. I am currently dissatisfied with this family, so as more exemplars come in, this set should eventually get recomputed. 

## Pneumoviridae
The cleaved forms of the fusion glycoprotein differ between ortho- and metapneumoviride.  All of the orthos, except murine and close relatives, have a p27 peptide that is a real protein. This necessitated the insertion of three taxon-level directories (ortho, meta, and murine orthos).  I kept the orginal all-pneumo alignment directory which has everything and has seprate subdirectories for the mature F proteins.   The pssm directories for the three taxa contain the pssms that I had originally built for all pneumos. This means that there are a few extra pssms that won't match and can be cleaned up on a rainy day.    






