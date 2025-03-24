# Current Protein Trees
The March 2025 alignments and trees have been computed for the *Bunyavirales*, *Filoviridae*, *Paramyxoviridae*, and *Pneumoviridae* using the annotation system.  *No* post-annotation curation has been performed on the alignments or trees. For each taxon, there are two subdirectories: "Good_Only" and "With_Poor".  The former includes only proteins from "Good" quality genomes, as determined by the annotation tool.  The latter, considers proteins from both "Good" and "Poor" quality genomes as long as there is no exception associated with a given feature in a poor quality genome (*e.g*, wrong length, ambiguous bases, or incorrect copy number).  For the *Bunyavirales*, there are two top-level directories.  One is the set of merged segments from BV-BRC, and the other is the set of segments unified by a single assembly accession downloaded from NCBI using their datasets tool. For BV-BRC, the segments of the *Bunyavirales* were linked using their full genome name string, which includes strain names and numbers.  All directories contain an html visualization of the tree, the Newick, aligned proteins in fasta format, and unaligned proteins.  Files are named based on annotation.<br> 

**Rules for making the trees** <br>
1.  Only segments >1Kb and <30Kb were considered.<br>
2.  No protein with an ambiguous X residue is allowed.<br>
3.  Alignments were computed with mafft.
4.  Trees were computed with FastTree.
5.  HTML trees were generated using FangFang's tool, svr_html (part of the SEED svr_script distrbution.
6.  HTML trees are colored by the family chosen for their annotation.

It is worth noting that occasionally very distant tips will be a higher BLAST match a taxon that isn't expected. This is more rare in the good quality genomes because they are full-length. I currently think the behavior is acceptable.



## Current Bunyavirales Protein Trees
### From BV-BRC, Only Proteins from Good Quality Genomes
[Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html)<br>

[Bunyavirales_small_nonstructural_protein_NSs_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Bunyavirales_small_nonstructural_protein_NSs_protein.html)<br>

[Fimoviridae_uncharacterized_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Fimoviridae_uncharacterized_protein.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html)<br>

[Movement_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Movement_protein.html)<br>

[Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Phenuiviridae_mature_nonstructural_78-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Phenuiviridae_mature_nonstructural_78-kD_protein.html)<br>

[Phenuiviridae_nonstructural_protein_5_NS5_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Phenuiviridae_nonstructural_protein_5_NS5_protein.html)<br>

[Phenuiviridae_Pc3_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Phenuiviridae_Pc3_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[RNA_silencing_suppressor_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/RNA_silencing_suppressor_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Signal_peptide.html)<br>

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br>





### From BV-BRC, Including Proteins from Poor Quality Genomes

[Arenaviridae_RING_finger_protein_Z_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Arenaviridae_RING_finger_protein_Z_protein.html)<br>

[Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html)<br>

[Bunyavirales_small_nonstructural_protein_NSs_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Bunyavirales_small_nonstructural_protein_NSs_protein.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P7.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Fimoviridae_RNA_silencing_suppressor_protein_P7.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P8.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Fimoviridae_RNA_silencing_suppressor_protein_P8.html)<br>

[Fimoviridae_uncharacterized_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Fimoviridae_uncharacterized_protein.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html)<br>

[Movement_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Movement_protein.html)<br>

[Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Nucleocapsid_protein_N_protein.html)<br>

[Phenuiviridae_mature_nonstructural_78-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Phenuiviridae_mature_nonstructural_78-kD_protein.html)<br>

[Phenuiviridae_nonstructural_protein_5_NS5_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Phenuiviridae_nonstructural_protein_5_NS5_protein.html)<br>

[Phenuiviridae_Pc3_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Phenuiviridae_Pc3_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[RNA_silencing_suppressor_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/RNA_silencing_suppressor_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Signal_peptide.html)<br>

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br>


### From NCBI, Only Proteins from Good Quality Genomes with a single Assembly Accession
[Arenaviridae_RING_finger_protein_Z_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Arenaviridae_RING_finger_protein_Z_protein.html)<br>

[Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html)<br>

[Bunyavirales_small_nonstructural_protein_NSs_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Bunyavirales_small_nonstructural_protein_NSs_protein.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P7.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_RNA_silencing_suppressor_protein_P7.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P8.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_RNA_silencing_suppressor_protein_P8.html)<br>

[Fimoviridae_uncharacterized_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_uncharacterized_protein.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html)<br>

[Movement_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Movement_protein.html)<br>

[Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Phenuiviridae_mature_nonstructural_78-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_mature_nonstructural_78-kD_protein.html)<br>

[Phenuiviridae_nonstructural_protein_5_NS5_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_nonstructural_protein_5_NS5_protein.html)<br>

[Phenuiviridae_Pc3_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_Pc3_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[RNA_silencing_suppressor_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/RNA_silencing_suppressor_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Signal_peptide.html)<br>

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br>


### From NCBI, Including Proteins from Poor Quality Genomes with a single Assembly Accession

[Arenaviridae_RING_finger_protein_Z_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Arenaviridae_RING_finger_protein_Z_protein.html)<br>

[Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Bunyavirales_mature_nonstructural_membrane_protein_NSm_protein.html)<br>

[Bunyavirales_small_nonstructural_protein_NSs_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Bunyavirales_small_nonstructural_protein_NSs_protein.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P7.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_RNA_silencing_suppressor_protein_P7.html)<br>

[Fimoviridae_RNA_silencing_suppressor_protein_P8.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_RNA_silencing_suppressor_protein_P8.html)<br>

[Fimoviridae_uncharacterized_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Fimoviridae_uncharacterized_protein.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Mature_secreted_nonstructural_mucin-like_protein_MLD_protein.html)<br>

[Movement_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Movement_protein.html)<br>

[Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Nairoviridae_mature_secreted_nonstructural_38-kD_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Phenuiviridae_mature_nonstructural_78-kD_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_mature_nonstructural_78-kD_protein.html)<br>

[Phenuiviridae_nonstructural_protein_5_NS5_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_nonstructural_protein_5_NS5_protein.html)<br>

[Phenuiviridae_Pc3_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Phenuiviridae_Pc3_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[RNA_silencing_suppressor_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/RNA_silencing_suppressor_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Signal_peptide.html)<br>

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br>


## Current Filoviridae Protein Trees

### Only Proteins from Good Quality Genomes
[Filoviridae_matrix_protein_protein_VP40.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Filoviridae_matrix_protein_protein_VP40.html)<br>

[Filoviridae_matrix_protein_VP40.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Filoviridae_matrix_protein_VP40.html)<br>

[Filoviridae_nucleocapsid_maturation_factor_VP24.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Filoviridae_nucleocapsid_maturation_factor_VP24.html)<br>

[Filoviridae_RNA-dependent_RNA_polymerase_cofactor_protein_VP35.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Filoviridae_RNA-dependent_RNA_polymerase_cofactor_protein_VP35.html)<br>

[Filoviridae_transcription_activator_protein_VP30.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Filoviridae_transcription_activator_protein_VP30.html)<br>

[Mature_delta_peptide_viroporin_and_enterotoxin.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Mature_delta_peptide_viroporin_and_enterotoxin.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_soluble_secreted_glycoprotein_sGP.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Mature_soluble_secreted_glycoprotein_sGP.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Signal_peptide.html)<br>

[Small_soluble_glycoprotein_ssGP_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Small_soluble_glycoprotein_ssGP_protein.html)<br>

[Soluble_secreted_glycoprotein_sGP_precursor.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Soluble_secreted_glycoprotein_sGP_precursor.html)<br>

### Including Proteins from Poor Quality Genomes
[Filoviridae_matrix_protein_VP40.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Filoviridae_matrix_protein_VP40.html)<br>

[Filoviridae_nucleocapsid_maturation_factor_VP24.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Filoviridae_nucleocapsid_maturation_factor_VP24.html)<br>

[Filoviridae_RNA-dependent_RNA_polymerase_cofactor_protein_VP35.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Filoviridae_RNA-dependent_RNA_polymerase_cofactor_protein_VP35.html)<br>

[Filoviridae_transcription_activator_protein_VP30.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Filoviridae_transcription_activator_protein_VP30.html)<br>

[Mature_delta_peptide_viroporin_and_enterotoxin.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Mature_delta_peptide_viroporin_and_enterotoxin.html)<br>

[Mature_envelope_glycoprotein_Gc_Gc_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Mature_envelope_glycoprotein_Gc_Gc_protein.html)<br>

[Mature_envelope_glycoprotein_Gn_Gn_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Mature_envelope_glycoprotein_Gn_Gn_protein.html)<br>

[Mature_soluble_secreted_glycoprotein_sGP.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Mature_soluble_secreted_glycoprotein_sGP.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Nucleocapsid_protein_N_protein.html)<br>

[Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Pre-glycoprotein_polyprotein_GP_complex_GPC_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Signal_peptide.html)<br>

[Small_soluble_glycoprotein_ssGP_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Small_soluble_glycoprotein_ssGP_protein.html)<br>

[Soluble_secreted_glycoprotein_sGP_precursor.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Soluble_secreted_glycoprotein_sGP_precursor.html)<br>