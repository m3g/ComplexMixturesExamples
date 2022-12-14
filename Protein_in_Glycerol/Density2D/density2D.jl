using Plots
using LaTeXStrings
using Formatting
using ComplexMixtures, PDBTools
script_dir = @__DIR__

# Default plot options (too look good with LaTeXStrings)
plot_font = "Computer Modern"
default(fontfamily=plot_font, linewidth=2, framestyle=:box, label=nothing)

# PDB file of the system simulated
pdb = readPDB("$script_dir/../Data/system.pdb")

# Load results of a ComplexMixtures run
R = load("$script_dir/../Data/results_glyc50.json")

# Inform which is the solute
protein = select(pdb, "protein")
solute = Selection(protein, nmols=1)

# Inform which is the solvent
glycerol = select(pdb, "resname GLYC")
solvent = Selection(glycerol, natomspermol=14)

# Plot a 2D map showing the contributions of some residues
residues = collect(eachresidue(protein))
irange = 70:110
rescontrib = zeros(length(R.mddf), length(residues))
for (ires, residue) in pairs(residues)
  rescontrib[:, ires] .= contrib(solute, R.solute_atom, residue)
end

# Plot only for distances within 1.5 and 3.5:
idmin = findfirst(d -> d > 1.5, R.d)
idmax = findfirst(d -> d > 3.5, R.d)

# Obtain pretty labels for the residues in the x-axis
labels = PDBTools.oneletter.(resname.(residues)) .* format.(resnum.(residues))

contourf(irange, R.d[idmin:idmax], rescontrib[idmin:idmax, irange],
  color=cgrad(:tempo), linewidth=0.1, linecolor=:black,
  colorbar=:none, levels=5,
  xlabel="Residue", ylabel=L"r/\AA",
  xticks=(irange, labels[irange]), xrotation=60,
  xtickfont=font(6, plot_font),
  size=(500, 280))

savefig("$script_dir/density2D.pdf")


