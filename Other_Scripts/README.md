## list_annos_from_pssms.pl 
This is a simple perl script that is dumps that reads the pssm directory and dumps the title fields.`Current.PSSM.annos.tab` is just the output taken on whatever date is posted on the file.<br>

## ali_to_pssm.pl
This is a janky script that will build a pssm from an alignment.  It has an option -t to add a title to the top of the file.  Note that the system is currently reading the title field to assign the annotation.  It requires by-hand editing of the title field to  sure you don't end up with 2 titles and it wants to append a "1" to the beginning of the title field.   I may eventually go back and clean this up.  It requires `BlastInterface.pm`<br>

##fasta-cluster-pssm.pl 
This perl script will create automatic clusters, alignments, and pssms for a set of proteins that are deemed to be isofunctional homologs.  PSSMs and alignments must be hand curated, but it reduces the workload. Here is an example usage:<br>
`fasta-cluster-pssm.pl -a "Pre-glycoprotein polyprotein GP complex, GPC" -p "Peribunyaviridae.GPC" -tmp "Peribunyaviridae-GPC"  <Peribunyaviridae.GPC.unique.len.aa`<br>
