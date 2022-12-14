using PDBTools
using ComplexMixtures
script_dir = @__DIR__

# The complete trajectory file can be downloaded from: 
# https://drive.google.com/file/d/1BuXJ8AjBeduMSD2CkDJLDNxAAD2QNNg6/view?usp=sharing
full_traj = isfile("../Test/trajectories/traj_Glyc.dcd")
if full_traj
    trajectory_file = "../Test/trajectories/traj_Glyc.dcd"
else
    println("WARNING: will execute calculations with a small trajectory sample.")
    trajectory_file = "../Test/trajectories/traj_Glyc_sample.dcd"
end

# Load system PDB file
system = readPDB("./simulation/equilibrated.pdb")

# Select the atoms corresponding to glycerol and water (using PDBTools)
glyc = select(system,"resname GLLM")
water = select(system,"water")

# Compute Glycerol-Glycerol auto correlation mddf 
solute = Selection(glyc,natomspermol=14)
traj = Trajectory(trajectory_file,solute) # solute and solvent are the same

# We define a large solute domain (large dbulk) to obtain a good convergence
# for the KB integral. The mddf converges at much shorter distances.   
opt = Options(dbulk=20.0)
mddf_glyc = mddf(traj,opt)

# Save results for later analysis
if full_traj
    save(mddf_glyc,"./results/mddf_glyc.json")
end

# Compute water-glycerol mddf
solvent = Selection(water,natomspermol=3)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=20.)
mddf_water_glyc = mddf(traj,opt)

# Save results for later analysis
if full_traj
    save(mddf_glyc,"./results/mddf_glyc.json")
end
