# Viral Annotation
This repo contains code and data for improving viral annotation.  It currently covers the members of the *Paramyxoviridae* and the *Arenaviridae*. The overall goal is to create a low-tech solution for calling viral proteins across entire viral families, and to cover cases where we do not have bespoke species-specific annotations.<br>

This program is not intended to be used as a *de novo* protein or ORF discovery tool.  It is designed to find proteins that we already know to exist.  I will explain more about how it works below.



## Dependencies

Unless otherwise stated, the programs described in this repo are written and tested in in perl (v5.38.0).

The script(s) have the following dependencies:<br>
External CPAN perl modules:

`JSON::XS`<br>
`File::Slurp`<br>

It also uses `gjoseqlib.pm` which is perl module that was written by Gary Olsen at the University of Illinois.  It is used for sequence manipulation.  You can get the latest version of this module by downloading the BV-BRC command line interface application. I have added a copy to this repo for convenience fromthe PATRIC Command Line Interface version 1.039 (the manipulation scripts are very stable, so I don't think anything I'm using has changed in a long time).  <br>

The program(s) run the blast suite of tools from NCBI.  The current version requires:<br>
`blastn: 2.13.0+`
`tblastn: 2.13.0+`<br>

I have not tried it on other versions of BLAST<br>

Although it should not matter, this has been developed using Bob's internal environment by sourcing:<br>
`/vol/patric3/cli/ubuntu-cli/user-env.sh`<br>

## How it works

The code is currently designed to work on *Paramyxoviridae* and *Arenaviridae* viruses, although I plan to add more.  As depicted in the image below, it first performs a blastn against a small set of representative genomes for each genus.  Then it sorts the results by bit score and chooses the best match.<br><br>
For each genus, there is a directory of PSSMs corresponding to each known protein for that genus. The PSSMs are derived from a set of hand curated alignments. In the next step, it cycles through each directory of PSSMs (there may be more than one PSSM per protein), choosing the best tblastn match per pssm. <br>

Note that it assumes your genome will have the same set of proteins as the nearest match. This is why it is not intended to be used as a discovery tool.  In the even that a new protein is found, a new PSSM must be added to the PSSM directory.  <br><br>

![Anno-Strategy](https://github.com/jimdavis1/Viral_Annotation/assets/7661533/0d6a3a44-47af-40bf-852d-5ddda250ad94)

<br><br>Finally it performs any special rules on the proteins/ORFs.  These rules are currently encoded in a JSON file called `Viral_PSSM.json`. The following is a description fo the current JSON strucutre. <br>

The perl code reads the JSON file into a hashref, which has this general structure: `options->{virus}->{protein}->{option} = value`.  So, for the F protein shown below, the options hash will have a bit score cutoff of 100 (fairly relaxed for a pssm), and a coverage cutoff of 65%. Upstream extension is turned on. This functionality extends the ORF upstream to find the closest Met (assuming it doesn't start with an AUG already). Downstream extension is also turned on.  This searches for a stop codon in-frame after the PSSM match. These can be turned off editing the JSON file and setting their values to zero.<br><br>

```
 "Metaavula": {
    "F": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1
    },
```


Another simple rule, which is not shown in the example is `start_to_met`: when this is set to 1, the first amino acid matching the PSSM is manually converted to Met.<br>

There are several cases where the rules are more complex.  For example, In the *Paramyxoviridae*, there is a phenomenon called RNA editing, where additional nucleotides can be inserted into transcripts of the phosphoprotein, which cause the RNA-polymerase to jump frames, and translation is continued in a new frame.  This means that two separate blast matches are required to capture the resulting protein. The functionality is encoded using the parameters `paramyxo_join`, `join_partner`, `new_anno`, and optionally `paramyxo_insert`. `paramyxo_join` is set to 1 if it is the first blast match of the pair and 2 if it is the second blast match of the pair. If `paramyxo_join` is set to 1, then the parameter `join_partner` must be used to denote the other pssm match that it is paired with.  In the case below, P is joined with V-ORF2 and W-ORF2 to to merge their ORFs and make a full-length V protein and full-length W protein, respectively.  On the second join partner, there is a parameter called `new_anno`, which carries the annotation for the newly merged protein (the annotations are typically stored in the pssm title field, but this could be changed). `paramyxo_insert` is a parameter that I have been tinkering with to correct the transcirpt, so that we end up with the right amino acid sequence after the merge.  Notice also that the `upstream_ext` parameter is set to zero in the case of V-ORF2.  In this case we do not want it searching up stream for a Met codon.  *Note that this phospho protein region is still a work in progress. Currently, I have tested this extensively and I have the gene boundaries debuged. In all cases that I am aware of it provides the merged amino acid product from the correct frame.  However, I still have several instances where the amino acid at the editing site is incorrect.  I also need to split the Respiroviruses*

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

The next and final set of rules relate to calling features that are not based on PSSMs.  In this case, if you know the location of a feature, and it can be established in reference to the start or stop position of any of the pssm-based matches, then it can be called in the genome.  Here is an example.  In this case for the *Arenaviridae* there is a feature called the large segement stemloop.  It ocurrs between the stop position of Z protein and the stop position of L protein. First, for L and Z we hadd a new key called `non_pssm_partner` which tells the program that one of the coordinates of that protein are important for calling the "Large Segment Stemloop."  Then in the object for the Large Segment Stemloop there are several new fields.  `anno` provides the annotation for the feature.  `translate` is a Boolean for whether you want the feature translated or not. `min_len` and `max_len` dictate the size cutoffs for calling the feature.  Next, `begin` tells it which pssm is the beginning coordinate with `begin_pssm`, `begin_pssm_loc` in this case is the "STOP" position of "Z".  Finally, `begin_offset` is 1 indicating that you want to start on the nucleotide after the stop.  `End` is essentially the mirror image.  If either of these had started or stoped on the "START"  position then you would simply say "START" for end or begin_pssm_loc.

```
  "Arenaviridae": {
    "L": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1,
      "non_pssm_partner": ["Large Segment Stemloop"]
    },
    "Z": {
      "bit_cutoff": 100,
      "coverage_cutoff": 0.65,
      "upstream_ext": 1,
      "downstream_ext": 1,
      "non_pssm_partner": ["Large Segment Stemloop"]
    },
    "Large Segment Stemloop": {
      "anno": "Large Segment Intergenic Stem Loop Region",
      "translate": 0,
      "min_len": 50,
      "max_len": 200,
      "begin": {
        "begin_pssm": "Z",
        "begin_pssm_loc": "STOP",
        "begin_offset": 1
      },
      "end": {
        "end_pssm": "L",
        "end_pssm_loc": "STOP",
        "end_offset": 1
      }
```




## Repo Contents




