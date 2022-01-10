# CABSdock protein-peptide global docking evaluation
CABS-dock standalone is a multiplatform Python package for protein–peptide docking with backbone flexibility. This python script is used to measure the IL-RMSD (RMSD of the ligand in the interface) proposed in Journal of Chemical Theory and Computation 2020 16 (6), 3959-3969, and the proportion of interface residues that was predicted correctly.

IL_RMSD was calculated on the basis of the backbone atoms of the peptide residues within 10 Å from the protein after the optimal superimposition of the protein residues within 10 Å from the peptide.
