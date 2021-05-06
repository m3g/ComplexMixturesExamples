using PDBTools
using ComplexMixtures
const CM = ComplexMixtures

# Load PDB file of the system
atoms = readPDB("../Data/system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms,"protein")
glyc = select(atoms,"resname GLYC")

# Setup solute and solvent structures
solute = CM.Selection(protein,nmols=1)
solvent = CM.Selection(glyc,natomspermol=14)

# Setup the Trajectory structure
trajectory = CM.Trajectory("../Data/glyc50.dcd",solute,solvent)

# Run the calculation and get results
results = CM.mddf(trajectory)

# Save the reults to recover them later if required
CM.save(results,"./glyc50.json")

# We will now load the "true" results file, computed from the 
# full trajectory. 
results = CM.load("../Data/results_glyc50.json")

# Produce plots
using Plots, Plots.PlotMeasures, LaTeXStrings
default(fontfamily="Computer Modern",
        linewidth=2, framestyle=:box, label=nothing, grid=false)

# The complete MDDF and the Kirkwood-Buff Integral
plot(layout=(1,2))
plot!(results.d,results.mddf,
      xlabel=L"r/\AA",ylabel="mddf",subplot=1)
hline!([1],linestyle=:dash,linecolor=:gray,subplot=1)
plot!(results.d,results.kb/1000, #to L/mol
      xlabel=L"r/\AA",ylabel=L"G_{us}/\mathrm{L~mol^{-1}}",subplot=2)
plot!(size=(800,300),margin=4mm)
savefig("./mddf.pdf")

# Atomic contributions to the MDDF
hydroxils = ["O1","O2","O3","H1","H2","H3"]
hydr_contrib = CM.contrib(solvent,results.solvent_atom,hydroxils)
aliphatic = ["C1","C2","HA","HB","HC","HD"]
aliph_contrib = CM.contrib(solvent,results.solvent_atom,aliphatic)

plot(results.d,results.mddf,xlabel=L"r/\AA",ylabel="mddf",size=(600,400))
plot!(results.d,hydr_contrib,label="Hydroxils")
plot!(results.d,aliph_contrib,label="Aliphatic chain")
hline!([1],linestyle=:dash,linecolor=:gray)
savefig("./mddf_atom_contrib.pdf")









