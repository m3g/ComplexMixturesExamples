using ComplexMixtures
using PDBTools

# Load pdb file of the system
system = readPDB("./simulation/equilibrated.pdb")

# The trajectory file can be downloaded from:
# https://drive.google.com/file/d/12TT5tblkFp1NtFOAQgjjGhmnYaXA8vQi/view?usp=sharing
trajectory_file = "../trajectories/traj_POPC.dcd"

# Select the atoms corresponding to glycerol and water
popc = select(system,"resname POPC")
water = select(system,"water")
ethanol = select(system,"resname ETOH")

# Set the complete membrane as the solute. We use nmols=1 here such
# that the membrane is considered a single solute in the calculation. 
solute = Selection(popc,nmols=1)

# Compute water-POPC distribution and KB integral 
solvent = Selection(water,natomspermol=3)
traj = Trajectory(trajectory_file,solute,solvent)

# We want to get reasonably converged KB integrals, which usually
# require large solute domains. Distribution functions converge 
# rapidly (~10Angs or less), on the other side.
opt = ComplexMixtures.Options(dbulk=20.)

# Compute the MDDF
mddf_water_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_water_POPC,"./results/mddf_water_POPC.json")

# Compute water-POPC distribution and KB integral 
solvent = Selection(ethanol,natomspermol=9)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=20.)
mddf_ethanol_POPC = mddf(traj,opt)

# Save results for later use
save(mddf_ethanol_POPC,"./results/mddf_ethanol_POPC.json")

