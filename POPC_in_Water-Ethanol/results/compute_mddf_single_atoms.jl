using ComplexMixtures
using PDBTools

# Load pdb file of the system
system = readPDB("./simulation/equilibrated.pdb")

# The trajectory file can be downloaded from:
# https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
trajectory_file = "../trajectories/traj_POPC.dcd"

# Select the atoms corresponding to glycerol and water
popc = select(system,"resname POPC")
O_water = select(system,"water and element O")
O_ethanol = select(system,"resname ETOH and element O")

# Set the complete membrane as the solute. We use nmols=1 here such
# that the membrane is considered a single solute in the calculation. 
solute = Selection(popc,nmols=1)

# Compute water-POPC distribution and KB integral 
solvent = Selection(O_water,natomspermol=1)
traj = Trajectory(trajectory_file,solute,solvent)

# We want to get reasonably converged KB integrals, which usually
# require large solute domains. Distribution functions converge 
# rapidly (~10Angs or less), on the other side.
opt = ComplexMixtures.Options(dbulk=20.,n_random_samples=2)

# Compute the MDDF
mddf_O_water_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_O_water_POPC,"./results/mddf_O_water_POPC.json")

# Compute water-POPC distribution and KB integral 
solvent = Selection(O_ethanol,natomspermol=1)
traj = Trajectory(trajectory_file,solute,solvent)
mddf_O_ethanol_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_O_ethanol_POPC,"./results/mddf_O_ethanol_POPC.json")

C2_ethanol = select(system,"resname ETOH and name C2")
solvent = Selection(C2_ethanol,natomspermol=1)
traj = Trajectory(trajectory_file,solute,solvent)
mddf_C2_ethanol_POPC = mddf(traj,opt)
save(mddf_C2_ethanol_POPC,"./results/mddf_C2_ethanol_POPC.json")





