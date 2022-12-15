using PDBTools
using ComplexMixtures
script_dir = @__DIR__

# PDB file of the system simulated
pdb = readPDB("$script_dir/../Data/system.pdb")

# Load results of a ComplexMixtures run
R = load("$script_dir/../Data/results_glyc50.json")

# Inform which is the solute
protein = select(pdb, "protein")
solute = Selection(protein, nmols=1)

# Compute the 3D density grid and output it to the PDB file
# here we use dmax=3.5 such that the the output file is not too large
grid = grid3D(
    solute=solute,
    solute_atoms=protein,
    mddf_result=R,
    output_file="$script_dir/grid.pdb",
    dmin=1.5,
    dmax=3.5
)





