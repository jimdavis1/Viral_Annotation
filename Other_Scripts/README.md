## list_annos_from_pssms.pl 
This is a simple perl script that passes over the pssm directory and dumps the title fields. `Current.PSSM.annos.tab` is just the output.<br>

## ali_to_pssm.pl
This is a janky script that will build a pssm from an alignment.  It has an option -t to add a title to the top of the file.  Note that the system is currently reading the title field to assign the annotation.  It requires by-hand editing of the title field to  sure you don't end up with 2 titles and it wants to append a "1" to the beginning of the title field.   I may eventually go back and clean this up.  It requires `BlastInterface.pm`<br>

## fasta-cluster-pssm.pl 
This perl script will create automatic clusters, alignments, and pssms for a set of proteins that are deemed to be isofunctional homologs.  PSSMs and alignments must be hand curated, but it reduces the workload.  It runs system calls to `mmseqs` and `mafft`. 
 It also uses BlastInterface. Here is an example usage:<br>
`fasta-cluster-pssm.pl -a "Pre-glycoprotein polyprotein GP complex, GPC" -p "Peribunyaviridae.GPC" -tmp "Peribunyaviridae-GPC"  <Peribunyaviridae.GPC.unique.len.aa`<br>

## viral_genome_quality.pl
This perl script assesses genome quality based on the presence or absence of key proteins in the genome by checking their copy number and lengths. It also calls the contigs (when possible) based on the presence of these key proteins, and evaluates the contigs based on their copy number, lengths, and fraction of ambiguous nucleotides.  It works by reading in a GTO of an annotated virus.  The genome must first be annotated by `annotate_by_viral_pssm-GTO.pl`, becuase it is specifically designed to only look for our annotatons.  Specifically, in each key taxonomic group it is looking for the following:<br>
1. Paramyxoviridae:  F, L, M, and N.  (The P protein is obviously essential, but because of the weirdness with polymerase stuttering, it is ignored).
2. Arenaviridae:  GPC, L, N, and Z (mat peptides are not currently considered for simplicity).
3. Hantaviridae, Nairoviridae, Peribunyaviridae: GPC, L, N.   In the Nairos, according to ICTV, there are a small number of metagenomically characterized tick-associated viruses that lack a M segment (and a known GPC).  It is unknown if the GPC was lost, or if there is some other protein doing the job of the GPC peptides.  Either way, if a 2-segment Nairovirus is encountered, it will be considered incomplete. 
4. Phenuiviridae: L and N.  Note: The plant viruses of this group do not have a GPC. According to ICTV they can have 2-8 segments. Becuase of this, the program only calls the "L-encoding RNA Segment" and the "N-encoding RNA Segment".<br>
<br>
Usage<br> `viral_genome_quality.pl -i Input.gto -o Output.gto -p Prefix`.<br>
The prefix is for the `Prefix.feature_quality` and `Prefix.contig_quality` output files.  STDOUT is `GenomeID\t"Good" or "Poor"`.
<br> 
