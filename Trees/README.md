# Current Protein Trees
The following trees were generated **without curation** from reasonably high quality genomes in BVBRC. The intention is to work toward having a system that can self-generate and edit the alignments and trees.<br> 

**Rules for making the trees** <br>
1.  An attempt was made to merge segments for the Bunyavirales based on genome name. In order to be annotated, a genome had to have at least 2 segments >= 1KB. For the Paramyxos, the contig length had to be >= 800nts.<br>
2.  No protein with an ambiguous X residue is allowed in the tree.<br>
3.  Alignments were automatically curated to cull sequences that were A) too short or B) had extremely aberant occupancy patterns. For (B) occupancy is computed by first determining the fraction occupancy for each column in the alignment.  These occupancy fractions are summed for each occupied column, for each sequence.  Any sequence >4 standard deviations away from the occupancy of all sequences or with an average within-sequence occupancy >4 standard deviations from the mean were eliminated.  These are listed in the .culled.fa files in this repo.  It turns out to be a fairly low bar.  
4.  Finally branch lengths in the tree exceeding 1.25x the standard deviation of the average branch lengths for the entire tree were pruned. <br>

All alginments, trees, and HTML files are available in subsequent directories. HTML files are colored by the family that they have been assigned based on the NCBI taxonomy from GenBank.


## Current Bunyavirales Protein Trees

[Bunya.78-kD_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.78-kD_Protein.html)<br>

[Bunya.Gc_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.Gc_Protein.html)<br>

[Bunya.Gn_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.Gn_Protein.html)<br>

[Bunya.GP38_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.GP38_Protein.html)<br>

[Bunya.GPC_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.GPC_Protein.html)<br>

[Bunya.L_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.L_Protein.html)<br>

[Bunya.MLD_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.MLD_Protein.html)<br>

[Bunya.Movement_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.Movement_Protein.html)<br>

[Bunya.N_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.N_Protein.html)<br>

[Bunya.NSm_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.NSm_Protein.html)<br>

[Bunya.NSs_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.NSs_Protein.html)<br>

[Bunya.SSP.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.SSP.html)<br>

[Bunya.Z_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Bunyavirales/Bunya.Z_Protein.html)<br>



## Current Paramyxoviridae Protein Trees

[Paramyxo.C-prime_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.C-prime_Protein.html)<br>

[Paramyxo.C_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.C_Protein.html)<br>

[Paramyxo.D_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.D_Protein.html)<br>

[Paramyxo.F_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.F_Protein.html)<br>

[Paramyxo.G_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.G_Protein.html)<br>

[Paramyxo.HN_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.HN_Protein.html)<br>

[Paramyxo.H_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.H_Protein.html)<br>

[Paramyxo.I_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.I_Protein.html)<br>

[Paramyxo.L_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.L_Protein.html)<br>

[Paramyxo.M_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.M_Protein.html)<br>

[Paramyxo.N_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.N_Protein.html)<br>

[Paramyxo.ORF1.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.ORF1.html)<br>

[Paramyxo.P_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.P_Protein.html)<br>

[Paramyxo.SH_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.SH_Protein.html)<br>

[Paramyxo.TM_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.TM_Protein.html)<br>

[Paramyxo.U_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.U_Protein.html)<br>

[Paramyxo.V_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.V_Protein.html)<br>

[Paramyxo.W_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.W_Protein.html)<br>

[Paramyxo.Y1_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.Y1_Protein.html)<br>

[Paramyxo.Y2_Protein.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/Paramyxoviridae/Paramyxo.Y2_Protein.html)<br>

## Current Whole Genome Sequence Trees
The Paramyxos are whole genome DNA trees, the Bunyavirales are a concatenation of gene trees for L, N, and GPC.  Note that the alignments are not hosted here since they are too big.

[Paramyxo.clean.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/WGS-Trees/Paramyxo.clean.html)<br>
[Bunya.L-N-GPC.clean.html](https://htmlpreview.github.io/?https://github.com/jimdavis1/Viral_Annotation/blob/main/Trees/June_2024/WGS-Trees/Bunya.L-N-GPC.clean.html)<br>


