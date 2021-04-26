# 2D residue contribution density map

In this example we compute the density map of Glycerol in the vicinity of a set of residues of a protein, from the minimum-distance distribution function. 

The MDDF can be decomposed in the contributions of each atom of the solute or of the solvent. Here, we sum up te contributions of all the atoms of each residue of the solute, which is a protein, and plot a density map with the final information. The output figure obtained is:

<center>
<img src="./density.png">
</center>

## How to run this example:

1. Download and install [Julia](https://julialang.org).

2. Install all required packages. Within Julia, do:
```julia
julia> ] add ComplexMixtures, PDBTools, Plots, LaTeXStrings, Formatting
```

3. Get all files: 
```bash
git clone https://github.com/m3g/ComplexMixturesExamples
```

4. Run the example:
```bash
cd ComplexMixturesExamples/Density
julia density.jl
```

## Detailed explanation of the example:

Here, we use the `contourf` function of the `Plots` package of Julia. A detailed explanation of the input file `density.jl` is provide below: 

### Loading packages that will be used:


```julia
using Plots
using LaTeXStrings
using Formatting
using ComplexMixtures, PDBTools
const CM = ComplexMixtures
```

### Some default options so the plot looks nice
```julia
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing)
```

### Read the PDB file (using PDBTools)
```julia
pdb = readPDB("./system.pdb")
```

### Load results of the ComplexMixtures run
```julia
R = CM.load("./results_glyc50.json")  
```

### Define which are the solute molecules (the protein)
```julia
protein = select(pdb,"protein")
solute = CM.Selection(protein,nmols=1)
```

### Define which are the solvent molecules (Glycerol here)
```julia
glycerol = select(pdb,"resname GLYC")
solvent = CM.Selection(glycerol,natomspermol=14)
```

### Retrive the resiude contribution data

Collect which are the protein residues 
```julia
residues = collect(eachresidue(protein))
```

Set a matrix that will store the results, with a number of lines corresponding to the length of the MDDF histogram, and with a number of columns corresponding to the number of residues:
```julia
rescontrib = zeros(length(R.mddf),length(residues))
```

Now, collect the contribution of each residue as a column of the above matrix. The notation `pairs(residues)` returns tuples containg the index `ires` and the corresponding residue. The `.=` symbol sets each element of the corresponding column of the  `rescontrib` matrix to the output of `CM.contrib` (by broadcasting).  
```julia
for (ires,residue) in pairs(residues)
  rescontrib[:,ires] .= CM.contrib(solute,R.solute_atom,residue)
end
```

### Plot only for distances within 1.5 and 3.5:

Here, we will plot only the contributions from residue `70` to residue `110`, and from distances ranging from `1.5` to `3.5` which is where most of the action occurs:
```julia
irange=70:110
idmin = findfirst( d -> d > 1.5, R.d)
idmax = findfirst( d -> d > 3.5, R.d)
```

To obtain pretty labels for the residues in the x-axis, we retrive the one-letter residue names and concatenate them with the residue number converted to strings:

```julia
labels = PDBTools.oneletter.(resname.(residues)).*format.(resnum.(residues))
```

And, finally, we produce the plot, with a series of options that make this particular contour plot look nice:

```julia
contourf(irange,R.d[idmin:idmax],rescontrib[idmin:idmax,irange],
         color=cgrad(:tempo),linewidth=0.1,linecolor=:black,
         colorbar=:none,levels=5,
         xlabel="Residue",ylabel=L"r/\AA",
         xticks=(irange,labels[irange]),xrotation=60,
         xtickfont=font(6,plot_font),
         size=(500,280))
```

The final figure is saved as a `pdf` file:
```julia
savefig("./density.pdf")
```




