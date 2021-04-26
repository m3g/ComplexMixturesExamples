# 2D residue contribution density map

In this example we compute the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

The MDDF can be decomposed in the contributions of each atom of the solute or of the solvent. Here, we sum up te contributions of all the atoms of each residue of the solute, which is a protein, and plot a density map with the final information. The output figure obtained is:

<center>
<img src="./density.png">
</center>

Here, we use the `contourf` function of the `Plots` package of Julia. A detailed explanation of the input file `density.jl` is provide below: 


