using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

# Plot defaults
plot_font = "Computer Modern"
default(fontfamily=plot_font,linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes()

system = readPDB("./equilibrated.pdb")

# The results will be loaded from previous computations. The original
# code for computing the distributions is commented.
#trajectory_file = "../trajectories/traj_Glyc.dcd"

# Select the atoms corresponding to glycerol and water
glyc = select(system,"resname GLLM")
water = select(system,"water")

# Compute glyc-glyc auto correlation 
solute = Selection(glyc,natomspermol=14)
#traj = Trajectory(trajectory_file,solute)

#opt = ComplexMixtures.Options(dbulk=20.)
#mddf_glyc = mddf(traj,opt)
#save(mddf_glyc,"mddf_glyc.json")

# Load previously computed data
mddf_glyc = load("./mddf_glyc.json")

# Plot the glycerol auto-correlation function and KB integral
plot(layout=(2,1))
x = mddf_glyc.d
y = movavg(mddf_glyc.mddf,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel="MDDF",subplot=1)
plot!(xlim=(0,20),subplot=1)

y = movavg(mddf_glyc.kb,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",subplot=2)
plot!(xlim=(0,20),subplot=2)
savefig("mddf_kb.png")

# Plot some group contributions to the MDDF
hydroxyls = ["O1","O2","O3","HO1","HO2","HO3"]
aliphatic = ["C1","C2","C3","H11","H12","H2","H31","H32"]

hydroxyl_contrib = contrib(solute,mddf_glyc.solvent_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_glyc.solvent_atom,aliphatic)

x = mddf_glyc.d
y = movavg(mddf_glyc.mddf,n=10).x
plot(x,y,label="Total")
y = movavg(hydroxyl_contrib,n=10).x
plot!(x,y,label="Hydroxyl contributions")
y = movavg(aliphatic_contrib,n=10).x
plot!(x,y,label="Aliphatic contributions")
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("mddf_glyc_groups.png")

# Interactions of glycerol with water
solvent = Selection(water,natomspermol=3)
#traj = Trajectory(trajectory_file,solute,solvent)
#mddf_water_glyc = mddf(traj,opt)
#save(mddf_water_glyc,"mddf_water_glyc.json")
mddf_water_glyc = load("./mddf_water_glyc.json")

x = mddf_water_glyc.d
y = movavg(mddf_water_glyc.mddf,n=10).x
plot(x,y,label="Total",subplot=1)

hydroxyl_contrib = contrib(solute,mddf_water_glyc.solute_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_water_glyc.solute_atom,aliphatic)

y = movavg(hydroxyl_contrib,n=10).x
plot!(x,y,label="Hydroxyl contributions")
y = movavg(aliphatic_contrib,n=10).x
plot!(x,y,label="Aliphatic contributions")
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
)
savefig("mddf_water_glyc.png")

# 2D plot of group contributions
groups = [
    select(system,"name O1 HO1"), # first hydroxyl
    select(system,"name C1 H11 H12"), # first CH2
    select(system,"name O2 HO2"), # second hydroxyl
    select(system,"name C2 H2"), # CH
    select(system,"name C3 H31 H32"), # second CH2
    select(system,"name O3 OH3") # third hydroxyl
] 
group_contrib = zeros(length(mddf_glyc.d),length(groups))
for (igroup, group) in pairs(groups)
    group_contrib[:,igroup] .= contrib(solute,mddf_glyc.solute_atom,group)
end

idmin = findfirst( d -> d > 1.5, mddf_glyc.d)
idmax = findfirst( d -> d > 3.0, mddf_glyc.d)
labels = [ "OH", L"\textrm{CH_2}", "OH", "CH", L"\textrm{CH_2}", "OH" ] 

contourf(
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xlabel="Group",ylabel=L"r/\AA",
    xticks=(1:length(groups),labels),xrotation=60,
)
savefig("map2D_glyc_glyc.png")





















