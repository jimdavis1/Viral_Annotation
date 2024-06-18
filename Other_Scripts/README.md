## list_annos_from_pssms.pl 
This is a simple perl script that is dumps that reads the pssm directory and dumps the title fields.`Current.PSSM.annos.tab` is just the output taken on whatever date is posted on the file.<br>

## ali_to_pssm.pl
This is a janky script that will build a pssm from an alignment.  It has an option -t to add a title to the top of the file.  It requires by-hand editing to make sure you don't end up with 2 titles and it wants to append a "1" to the beginning of the title field.   I may eventually go back and clean this up.  It requires `BlastInterface.pm`
