using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

# Plot defaults
plot_font = "Computer Modern"
default(fontfamily=plot_font,linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(); scalefontsizes(1.3)

# Load pdb file of the sytem
system = readPDB("./simulation/equilibrated.pdb")

# The results will be loaded from previous computations. The original
# code for computing the distributions is commented.
trajectory_file = "../trajectories/traj_POPC.dcd"

# Select the atoms corresponding to glycerol and water
popc = select(system,"resname POPC")
water = select(system,"water")
ethanol = select(system,"resname ETOH")

# Compute water-POPC distribution and KB integral 
solute = Selection(popc,nmols=1)
solvent = Selection(water,natomspermol=3)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=10.,lastframe=200)
mddf_water_POPC = mddf(traj,opt)
#save(mddf_water_POPC,"./results/mddf_water_POPC.json")

# Compute water-POPC distribution and KB integral 
solute = Selection(popc,nmols=1)
solvent = Selection(ethanol,natomspermol=9)
traj = Trajectory(trajectory_file,solute,solvent)
opt = ComplexMixtures.Options(dbulk=10.,lastframe=200)
mddf_ethanol_POPC = mddf(traj,opt)
#save(mddf_ethanol_POPC,"./results/mddf_ethanol_POPC.json")

# Load previously computed data
mddf_water_POPC = load("./results/mddf_water_POPC.json")
mddf_ethanol_POPC = load("./results/mddf_ethanol_POPC.json")

# Water-POPC MDDF and KB integral
plot(layout=(2,1))
x = mddf_water_POPC.d
y = movavg(mddf_water_POPC.mddf,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel="MDDF",subplot=1)
plot!(xlim=(0,10),subplot=1)

y = movavg(mddf_water_POPC.kb,n=10).x/1000
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~L~mol^{-1}}",subplot=2)
plot!(xlim=(0,10),subplot=2)

# Ethanol-POPC MDDF and KB integral
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel="MDDF",subplot=1)
plot!(xlim=(0,10),subplot=1)

y = movavg(mddf_ethanol_POPC.kb,n=10).x/1000
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~L~mol^{-1}}",subplot=2)
plot!(xlim=(0,10),subplot=2)

savefig("./results/mddf_kb.png")

# Contributions of the ethanol groups
groups = [
    (select(ethanol,"name O1 or name HO1"),"Hydroxyl"),
    (select(ethanol,"not name O1 and not name HO1"),"Aliphatic chain"),
]
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot(x,y,label="Total")
for group in groups
    group_contrib = contrib(solvent,mddf_ethanol_POPC.solvent_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2])
end
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./results/mddf_ethanol_groups.png")

# Contributions of POPC groups to ethanol-POPC distribution
groups = [
  (["N","C12","H12A","H12B","C13","H13A","H13B","H13C","C14",
    "H14A","H14B","H14C","C15","H15A","H15B","H15C"],"Amine"),
  ([ "C11","H11A","H11B","P","O13","O14","O12","O11",
     "C1","HA","HB" ],"Phosphate"),
  (["C2","HS","O21","C21","O22","C22","H2R","H2S","C3",
    "HX","HY","O31","C31","O32","C32","H2X","H2Y"],"Esthers"),
  (["C23","H3R","H3S","C24","H4R","H4S","C25","H5R","H5S",
    "C26","H6R","H6S","C27","H7R","H7S","C28","H8R","H8S",
    "C29","H91","C210","H101","C211","H11R","H11S","C212",
    "H12R","H12S","C213","H13R","H13S","C214","H14R","H14S",
    "C215","H15R","H15S","C216","H16R","H16S","C217","H17R",
    "H17S","C218","H18R","H18S","H18T","C33","H3X","H3Y",
    "C34","H4X","H4Y","C35","H5X","H5Y","C36","H6X","H6Y",
    "C37","H7X","H7Y","C38","H8X","H8Y","C39","H9X","H9Y",
    "C310","H10X","H10Y","C311","H11X","H11Y","C312","H12X",
    "H12Y","C313","H13X","H13Y","C314","H14X","H14Y","C315",
    "H15X","H15Y","C316","H16X","H16Y", "H16Z" ], "Aliphatic chain")
]
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot(x,y,label="Total")
for group in groups
    group_contrib = contrib(solute,mddf_ethanol_POPC.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2])
end
plot!(
    xlim=(1.3,5),
    ylim=(0,2),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./results/mddf_ethanol_groups.png")

# the same for POPC-water
x = mddf_water_POPC.d
y = movavg(mddf_water_POPC.mddf,n=10).x
plot(x,y,label="Total")
for group in groups
    group_contrib = contrib(solute,mddf_water_POPC.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2])
end
plot!(
    xlim=(1.3,5),
    ylim=(0,2),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./results/mddf_water_groups.png")
