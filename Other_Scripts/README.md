# Other Scripts
This directory contains a set of helper scripts that I use for a variety of tasks that are related to viral annotation generation and curation. They are mostly here for my own personal convenience, but you may find them helpful.  I offer no guarantees, they are usually in a state of flux.  Perhaps the most useful file in this directory is **current_annos.tab** which contains the exhaustive list of all annotations that have been generated to date with their associated taxa. 


## Some useful program descriptions
### list_annos.pl 
This is a simple perl script that passes over the json options file and returns every annotation for every family. `current.annos.tab` is the output.<br>

### list_pssms.pl 
This is a simple perl script that passes over the pssm directory and dumps the title fields so that you can see every PSSM for every family.  Note that the annotation filed is deprecated.  I used to use this for the annotation but this has been moved ot the json file.<br>

### ali_to_pssm.pl
This is a janky script that will build a pssm from an alignment.  It has an option -t to add a title to the top of the file.  Note that the system is currently reading the title field to assign the annotation.  It requires by-hand editing of the title field to  sure you don't end up with 2 titles and it wants to append a "1" to the beginning of the title field.   I may eventually go back and clean this up.  It requires `BlastInterface.pm`<br>

### fasta-cluster-pssm.pl 
This perl script will create automatic clusters, alignments, and pssms for a set of proteins that are deemed to be isofunctional homologs.  PSSMs and alignments must be hand curated, but it reduces the workload.  It runs system calls to `mmseqs` and `mafft`. 
 It also uses BlastInterface. Here is an example usage:<br>
`fasta-cluster-pssm.pl -a "Pre-glycoprotein polyprotein GP complex, GPC" -p "Peribunyaviridae.GPC" -tmp "Peribunyaviridae-GPC"  <Peribunyaviridae.GPC.unique.len.aa`<br>

###annotate_viral_taxon_ContigDir
This program will run the annotation pipleline over a directory of contigs.  It requires a metadata file that is formatted as [file_name \t genome_name \t taxon_id].<br>

###process_quality_dir.pl
This program provides some rough stats on a set of genomes in a Quality_GTO directory.<br>

###gto-dir-to-alis-and-trees.pl
This program will grind out the alignments and trees for a directory of GTOs.<br>








