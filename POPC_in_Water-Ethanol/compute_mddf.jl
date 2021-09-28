using ComplexMixtures
using PDBTools

# Load pdb file of the sytem
system = readPDB("./simulation/equilibrated.pdb")

# The trajectory file can be downloaded from:
# https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
trajectory_file = "../trajectories/traj_POPC.dcd"

# Select the atoms corresponding to glycerol and water
popc = select(system,"resname POPC")
water = select(system,"water")
ethanol = select(system,"resname ETOH")

# Set the complete membrane as the solute
solute = Selection(popc,nmols=1)

# Compute water-POPC distribution and KB integral 
solvent = Selection(water,natomspermol=3)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=10.,lastframe=200)
mddf_water_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_water_POPC,"./results/mddf_water_POPC.json")

# Compute water-POPC distribution and KB integral 
solvent = Selection(ethanol,natomspermol=9)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=10.,lastframe=200)
mddf_ethanol_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_ethanol_POPC,"./results/mddf_ethanol_POPC.json")