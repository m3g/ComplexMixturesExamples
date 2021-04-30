using Plots
using LaTeXStrings
using Formatting
using ComplexMixtures, PDBTools
const CM = ComplexMixtures

# Default plot options (too look good with LaTeXStrings)
plot_font = "Computer Modern"
default(fontfamily=plot_font,
        linewidth=2, framestyle=:box, label=nothing)

# PDB file of the system simulated
pdb = readPDB("../Data/system.pdb")

# Load results of a ComplexMixtures run
R = CM.load("../Data/results_glyc50.json")  

# Inform which is the solute
protein = select(pdb,"protein")
solute = CM.Selection(protein,nmols=1)

# Inform which is the solvent
glycerol = select(pdb,"resname GLYC")
solvent = CM.Selection(glycerol,natomspermol=14)

#include("./mddf3D.jl")

#mddf_setgrid3D(solute,protein,R,)



