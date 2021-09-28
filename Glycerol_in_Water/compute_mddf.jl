using ComplexMixtures
using PDBTools

# Load system PDB file
system = readPDB("./simulation/equilibrated.pdb")

# The complete trajectory file can be downloaded from: 
# https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing
trajectory_file = "../trajectories/traj_Glyc.dcd"

# Select the atoms corresponding to glycerol and water (using PDBTools)
glyc = select(system,"resname GLLM")
water = select(system,"water")

# Compute Glycerol-Glycerol auto correlation mddf 
solute = Selection(glyc,natomspermol=14)
traj = Trajectory(trajectory_file,solute) # solute and solvent are the same
opt = ComplexMixtures.Options(dbulk=20.)
mddf_glyc = mddf(traj,opt)

# Save results for later analysis
save(mddf_glyc,"./results/mddf_glyc.json")

# Compute water-glycerol mddf
solvent = Selection(water,natomspermol=3)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=20.)
mddf_water_glyc = mddf(traj,opt)

# Save results for later analysis
save(mddf_glyc,"./results/mddf_glyc.json")