# Viral Annotation
This repo contains code and data for improving viral annotation.  It currently covers the members of the Paramyxoviridae and the Arenaviridae. The overall goal is to create a low-tech solution for calling viral proteins across entire viral families, rather than only offering species-specific annotation.<br>

Unless otherwise stated, the programs described in this repo are written and tested in in perl (v5.38.0).

## Dependencies
The script(s) have the following dependencies:<br>
External CPAN perl modules:

`JSON::XS`<br>
`File::Slurp`<br>

It also uses the following perl module that was written by Gary Olsen at the University of Illinois.  It is used for sequence manipulation.  You can get this module by downloading the BV-BRC command line interface application.  
`gjoseqlib`<br>

The program(s) run the blast suite of tools from NCBI.  The current version runs:<br>
`blastn: 2.13.0+`
`tblastn: 2.13.0+`<br>

Although it should not matter, this has been developed using Bob's internal environment by sourcing:<br>
`/vol/patric3/cli/ubuntu-cli/user-env.sh`<br>

## How it works

The code is currently designed to work on Paramyxoviridae and Arenaviridae viruses.  As depicted in the image below, it first performs a blastn against a small set of representative genomes for each genus.  Then it sorts the results by bit score and chosing the best match.<br><br>
For each genus, there is a directory of PSSMs corresponding to each known protein for that genus. The PSSMs are derived from a set of hand curated alignments. In the next step, it cycles through each directory of PSSMs, choosing the best tblastn match per pssm. <br><br>

![Anno-Strategy](https://github.com/jimdavis1/Viral_Annotation/assets/7661533/0d6a3a44-47af-40bf-852d-5ddda250ad94)

<br><br>Finally it performs any special rules on the proteins/ORFs.  These rules are currently encoded in a JSON file. In the snippet of JSON shown below, the perl code reads this as a hashref which has this general structure: `options->{virus}->{protein}->{option} = value`.  So, for the F protein, the options hash will have a bit score cutoff of 100 (fairly relaxed for a pssm), and a coverage cutoff of 65%.  Upstream extension is turned on. This functionality extends the ORF upstream to find the closest Met (assuming it doesn't start with an AUG already). Downstream extension is also turned on.  This searches for a stop codon in-frame after the PSSM match. These can be turned off by setting them to zero for a given protein.<br><br>

Other simple rules, which are not shown include `start_to_met`: when this is turned on, the first matching amino acid is manually converted to Met. *I should probably update the code here to only enable start-to-met functionality for feasible start codons. In the paramyxos, the only ones that I have encountered so far are Thr.*<br><br>
In the Paramyxoviridae, there is a phenomenon called RNA editing, where additional nucleotides can be inserted into transcripts of the phosphoprotein, which cause the RNA-polymerase to jump frames, and translation is continued in a new frame.  This means that two separate blast matches are required to capture the resulting protein. The functionality is encoded using the parameters `paramyxo_join`, `join_partner`, `new_anno`, and optionally `paramyxo_insert`. `paramyxo_join` is set to 1 if it is the first blast match of the pair and 2 if it is the second blast match of the pair. If `paramyxo_join` is set to 1, then the parameter `join_partner` must be used to denote the other pssm match that it is paired with.  In the case below, P is joined with V-ORF2 and W-ORF2 to make full-length V protein and full-length W protein.  On the second join partner, there is a parameter called `new_anno`, which carries the annotation for the newly merged protein. `paramyxo_insert` is a parameter that I have been using to correct the transcirpt, so that we end up with the right amino acid sequence after the merge.  *This is a work in progress*. 

Here is a snippet to describe the overall structure:<br>
```
 "Metaavula": {
    "F": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1
    },
    "P": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1,
      "paramyxo_join": 1,
      "join_partner": ["V-ORF2", "W-ORF2"]
    },
    "V-ORF2": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 0,
      "downstream_ext": 1,
      "paramyxo_join": 2,
      "new_anno": "V Protein"
```




