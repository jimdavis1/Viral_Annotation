# Current Protein Trees
The March 2025 alignments and trees have been computed for the *Bunyavirales*, *Filoviridae*, *Paramyxoviridae*, and *Pneumoviridae* using the annotation system.<br>

### Whole Genome Trees
WGS trees have been computed for the high quality single-segmented viruses: <br>

[Filoviridae](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/WGS_Trees/Filo_WGS.html)<br>

[Pneumoviridae](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/WGS_Trees/Pneumo_WGS.html)<br>

[paramyxoviridae](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/WGS_Trees/Paramyxo_WGS.html)<br>
<br>

### Protein Trees
Each directory contains an html visualization of the tree, the Newick file, aligned proteins in fasta format, and unaligned proteins in fasta format.  Files are named based on annotation.  *No* post-annotation curation has been performed on the alignments or trees. For each taxon, there are two subdirectories: "Good_Only" and "With_Poor".  The former includes only proteins from "Good" quality genomes, as determined by the annotation tool.  The latter considers proteins from both "Good" and "Poor" quality genomes as long as there is no exception associated with a given protein in a poor quality genome (*e.g*, wrong length, ambiguous bases, or incorrect copy number).<br>  

The *Bunyavirales* are multi-segmented, and the viral community has a bad habit of submitting each segment separately to GenBank as if it were its own genome. In order to assess genome quality, a whole genome must be considered. We attempted to merge Bunyavirales segments to creat *potentially* complete genomes in two ways.  First, the segments in BV-BRC, were linked using the full genome name string, which includes the strain names and numbers.  We also downloaded every Bunyavirales genome NCBI where the segments were unified under a single assembly accession, using their datasets tool.  These genome sets were analyzed separately, although most of the NCBI genomes exist in the BV-BRC set.<br> 


**Rules for making the trees** <br>
1.  Only segments >1Kb and <30Kb were considered.<br>
2.  No protein with an ambiguous X residue is allowed.<br>
3.  Alignments were computed with mafft.
4.  Trees were computed with FastTree.
5.  HTML trees were generated using FangFang's tool, svr_html (part of the SEED svr_script distribution).
6.  HTML trees are colored by the family chosen for their annotation.

It is worth noting that occasionally very distant tips will be a higher BLAST match to a taxon that isn't expected, resulting in its misidentification by the annotation tool. This is more common in the poor quality genomes due to their poorer sequence quality and shorter lengths. The current BLAST-based taxon identification represents a trade-off between being able to accurately identify and analyze well-characterized genomes, versus being permissive enough to provide some annotations to more distant relatives.  Presently, I think the balance is acceptable.  Overall, this happens most frequently in the tips of the paramyxos and pneumos, which used to be classified together in the same family.


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

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br><br>





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

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_BVBRC/With_Poor/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br><br>




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

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br><br>




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

[Small_Nonstructural_Protein_NSs_NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Bunyavirales_NCBI/Good_Only/Small_Nonstructural_Protein_NSs_NSs_Protein.html)<br><br>





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

[Soluble_secreted_glycoprotein_sGP_precursor.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/Good_Only/Soluble_secreted_glycoprotein_sGP_precursor.html)<br><br>




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

[Soluble_secreted_glycoprotein_sGP_precursor.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Filoviridae/With_Poor/Soluble_secreted_glycoprotein_sGP_precursor.html)<br><br>




## Current Paramyxoviridae Protein Trees

### Only Proteins from Good Quality Genomes
[Fusion_glycoprotein_F_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Fusion_glycoprotein_F_protein.html)<br>

[Matrix_protein_M_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Matrix_protein_M_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Paramyxoviridae_C-prime_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_C-prime_protein.html)<br>

[Paramyxoviridae_C_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_C_protein.html)<br>

[Paramyxoviridae_I_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_I_protein.html)<br>

[Paramyxoviridae_small_hydorophobic_protein_SH_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_small_hydorophobic_protein_SH_protein.html)<br>

[Paramyxoviridae_transmembrane_protein_TM_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_transmembrane_protein_TM_protein.html)<br>

[Paramyxoviridae_uncharacterized_protein_U_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_uncharacterized_protein_U_protein.html)<br>

[Paramyxoviridae_V_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_V_protein.html)<br>

[Paramyxoviridae_W_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_W_protein.html)<br>

[Paramyxoviridae_Y1_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_Y1_protein.html)<br>

[Paramyxoviridae_Y2_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Paramyxoviridae_Y2_protein.html)<br>

[Phosphoprotein_P_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Phosphoprotein_P_protein.html)<br>

[Receptor_binding_protein_RBP.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/Receptor_binding_protein_RBP.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br><br>



### Including Proteins from Poor Quality Genomes
[Fusion_glycoprotein_F_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Fusion_glycoprotein_F_protein.html)<br>

[Matrix_protein_M_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Matrix_protein_M_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Nucleocapsid_protein_N_protein.html)<br>

[Paramyxoviridae_C-prime_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_C-prime_protein.html)<br>

[Paramyxoviridae_C_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_C_protein.html)<br>

[Paramyxoviridae_I_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_I_protein.html)<br>

[Paramyxoviridae_small_hydorophobic_protein_SH_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_small_hydorophobic_protein_SH_protein.html)<br>

[Paramyxoviridae_transmembrane_protein_TM_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_transmembrane_protein_TM_protein.html)<br>

[Paramyxoviridae_uncharacterized_protein_U_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_uncharacterized_protein_U_protein.html)<br>

[Paramyxoviridae_V_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_V_protein.html)<br>

[Paramyxoviridae_W_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_W_protein.html)<br>

[Paramyxoviridae_Y1_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_Y1_protein.html)<br>

[Paramyxoviridae_Y2_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Paramyxoviridae_Y2_protein.html)<br>

[Phosphoprotein_P_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Phosphoprotein_P_protein.html)<br>

[Receptor_binding_protein_RBP.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/Receptor_binding_protein_RBP.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Paramyxoviridae/With_Poor/RNA-dependent_RNA_polymerase_L_protein.html)<br>




## Current Pneumoviridae Protein Trees

[Fusion_glycoprotein_precursor_F0_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Fusion_glycoprotein_precursor_F0_protein.html)<br>

[Matrix_protein_M_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Matrix_protein_M_protein.html)<br>

[Mature_C-terminal_fusion_protein_subunit_F1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Mature_C-terminal_fusion_protein_subunit_F1.html)<br>

[Mature_N-terminal_fusion_protein_subunit_F2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Mature_N-terminal_fusion_protein_subunit_F2.html)<br>

[Mature_p27_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Mature_p27_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Nucleocapsid_protein_N_protein.html)<br>

[Phosphoprotein_P_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Phosphoprotein_P_protein.html)<br>

[Pneumoviridae_nonstructural_protein_NS1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Pneumoviridae_nonstructural_protein_NS1.html)<br>

[Pneumoviridae_nonstructural_protein_NS2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Pneumoviridae_nonstructural_protein_NS2.html)<br>

[Pneumoviridae_small_hydorophobic_protein_SH_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Pneumoviridae_small_hydorophobic_protein_SH_protein.html)<br>

[Receptor_binding_protein_G_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Receptor_binding_protein_G_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Signal_peptide.html)<br>

[Transcription_elongation_factor_M2-1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Transcription_elongation_factor_M2-1.html)<br>

[Transcription-replication_factor_M2-2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/Good_Only/Transcription-replication_factor_M2-2.html)<br><br>



### Including Proteins from Poor Quality Genomes
[Fusion_glycoprotein_precursor_F0_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Fusion_glycoprotein_precursor_F0_protein.html)<br>

[Hypothetical_P2_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Hypothetical_P2_protein.html)<br>

[Matrix_protein_M_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Matrix_protein_M_protein.html)<br>

[Mature_C-terminal_fusion_protein_subunit_F1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Mature_C-terminal_fusion_protein_subunit_F1.html)<br>

[Mature_N-terminal_fusion_protein_subunit_F2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Mature_N-terminal_fusion_protein_subunit_F2.html)<br>

[Mature_p27_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Mature_p27_protein.html)<br>

[Nucleocapsid_protein_N_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Nucleocapsid_protein_N_protein.html)<br>

[Phosphoprotein_P_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Phosphoprotein_P_protein.html)<br>

[Pneumoviridae_nonstructural_protein_NS1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Pneumoviridae_nonstructural_protein_NS1.html)<br>

[Pneumoviridae_nonstructural_protein_NS2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Pneumoviridae_nonstructural_protein_NS2.html)<br>

[Pneumoviridae_small_hydorophobic_protein_SH_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Pneumoviridae_small_hydorophobic_protein_SH_protein.html)<br>

[Receptor_binding_protein_G_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Receptor_binding_protein_G_protein.html)<br>

[RNA-dependent_RNA_polymerase_L_protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/RNA-dependent_RNA_polymerase_L_protein.html)<br>

[Signal_peptide.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Signal_peptide.html)<br>

[Transcription_elongation_factor_M2-1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Transcription_elongation_factor_M2-1.html)<br>

[Transcription-replication_factor_M2-2.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/March_2025/Pneumoviridae/With_Poor/Transcription-replication_factor_M2-2.html)<br><br>
