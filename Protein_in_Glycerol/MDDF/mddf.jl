using PDBTools
using ComplexMixtures
script_dir = @__DIR__

# The complete trajectory file can be downloaded from (3Gb):
# https://drive.google.com/file/d/14M30jDHRwUM77hzbDphgbu8mcWFBcQrX/view?usp=sharing
full_traj = isfile("$script_dir/../../test/trajectories/glyc50_complete.dcd")
if full_traj
    trajectory_file = "$script_dir/../../test/trajectories/glyc50_complete.dcd"
else
    trajectory_file = "$script_dir/../../test/trajectories/glyc50_sample.dcd"
    println("WARNING: short trajectory sample: $(normpath(trajectory_file))")
end

# Load PDB file of the system
atoms = readPDB("$script_dir/../Data/system.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname GLYC")

# Setup solute and solvent structures
solute = Selection(protein, nmols=1)
solvent = Selection(glyc, natomspermol=14)

# Setup the Trajectory structure
trajectory = Trajectory(trajectory_file, solute, solvent)

# Run the calculation and get results
results = mddf(trajectory)

# Save the reults to recover them later if required
if full_traj
    save(results, "$script_dir/../Data/results_glyc50.json")
end

# Produce plots
using Plots, Plots.PlotMeasures, LaTeXStrings
default(fontfamily="Computer Modern", linewidth=2, framestyle=:box, label=nothing, grid=false)

# The complete MDDF and the Kirkwood-Buff Integral
plot(layout=(1, 2))
plot!(results.d, results.mddf,
    xlabel=L"r/\AA", ylabel="mddf", subplot=1)
hline!([1], linestyle=:dash, linecolor=:gray, subplot=1)
plot!(results.d, results.kb / 1000, #to L/mol
    xlabel=L"r/\AA", ylabel=L"G_{us}/\mathrm{L~mol^{-1}}", subplot=2)
plot!(size=(800, 300), margin=4mm)
savefig("$script_dir/mddf.png")
savefig("$script_dir/mddf.pdf")

# Atomic contributions to the MDDF
hydroxyls = ["O1", "O2", "O3", "H1", "H2", "H3"]
aliphatic = ["C1", "C2", "HA", "HB", "HC", "HD"]
hydr_contrib = contributions(solvent, results.solvent_atom, hydroxyls)
aliph_contrib = contributions(solvent, results.solvent_atom, aliphatic)

plot(results.d, results.mddf, xlabel=L"r/\AA", ylabel="mddf", size=(600, 400))
plot!(results.d, hydr_contrib, label="Hydroxyls")
plot!(results.d, aliph_contrib, label="Aliphatic chain")
hline!([1], linestyle=:dash, linecolor=:gray)
savefig("$script_dir/mddf_atom_contrib.png")
savefig("$script_dir/mddf_atom_contrib.pdf")