using PDBTools
using ComplexMixtures
script_dir = @__DIR__

# The full trajectory file is available at: 
# https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing
full_traj = isfile("$script_dir/../Test/trajectories/simulation/traj_Polyacry.dcd")
if full_traj
    trajectory_file = "$script_dir/../Test/trajectories/simulation/traj_Polyacry.dcd"
else
    trajectory_file = "$script_dir/../Test/trajectories/simulation/traj_Polyacry_sample.dcd"
end

# Load a PDB file of the system
system = readPDB("$script_dir/simulation/equilibrated.pdb")

# Select the atoms corresponding DMF molecules
dmf = select(system,"resname DMF")

# Select the atoms corresponding to the Poly-acrylamide
acr = select(system,"resname FACR or resname ACR or resname LACR")

# Set the solute and the solvent selections for ComplexMixtures
solute = Selection(acr,nmols=1)
solvent = Selection(dmf,natomspermol=12)

# Set the trajectory structure
traj = Trajectory(trajectory_file,solute,solvent)

# Use a large dbulk distance for better KB convergence
opt = ComplexMixtures.Options(dbulk=20.)

# Compute the mddf and associated properties
mddf_dmf_acr = mddf(traj,opt)

# Save results to file for later use
if full_traj
    save(mddf_dmf_acr,"$script_dir/results/mddf_dmf_acr.json")
end









