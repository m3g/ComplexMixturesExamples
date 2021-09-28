using PDBTools
using ComplexMixtures

# Load a PDB file of the system
system = readPDB("./simulation/equilibrated.pdb")

# Trajectory file
# The trajectory file is available at: 
# https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing
trajectory_file = "../trajectories/traj_Polyacry.dcd"

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
save(mddf_dmf_acr,"./results/mddf_dmf_acr.json")








