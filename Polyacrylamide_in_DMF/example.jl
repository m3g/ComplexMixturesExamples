using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

# Plot defaults
plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    palette=:tab10
)
scalefontsizes(); scalefontsizes(1.3)

system = readPDB("./simulation/equilibrated.pdb")

# The results will be loaded from previous computations. The original
# code for computing the distributions is commented.
trajectory_file = "../trajectories/traj_Polyacry.dcd"

# Select the atoms corresponding to glycerol and water
dmf = select(system,"resname DMF")
acr = select(system,"resname FACR or resname ACR or resname LACR")

# Compute glyc-glyc auto correlation 
solute = Selection(acr,nmols=1)
solvent = Selection(dmf,natomspermol=12)
#traj = Trajectory(trajectory_file,solute,solvent)
#opt = ComplexMixtures.Options(dbulk=20.)
#mddf_dmf_acr = mddf(traj,opt)
#save(mddf_dmf_acr,"./results/mddf_dmf_acr.json")

# Load previously computed data
mddf_dmf_acr = load("./results/mddf_dmf_acr.json")

# Plot the MDDF of DMF relative to PolyACR and its corresponding KB integral
plot(layout=(2,1))
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel="MDDF",subplot=1)
plot!(xlim=(0,20),subplot=1)

y = movavg(mddf_dmf_acr.kb,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",subplot=2)
plot!(xlim=(0,20),subplot=2)
savefig("./results/mddf_kb.png")

# Plot DMF group contributions to the MDDF
groups = [ 
    (["O"],"O"),
    (["HA"],"HA"),
    (["N"],"N"),
    (["CC","CT","HC1","HC2","HC3","HT1","HT2","HT3"],"Methyl groups")
]
plot(layout=(2,1))
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,label="Total",subplot=1)
for group in groups
    group_contrib = contrib(solvent,mddf_dmf_acr.solvent_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=1)
end
plot!(
    xlim=(1,5),
    ylabel="MDDF",
    subplot=1
)

# Plot ACR group contributions to the MDDF
groups = [ 
    (["OE1","CD"],"CO"),
    (["NE2","HE22","HE21"],L"\textrm{NH_2}"),
    (["C","H2","H1","CA","HA","CF","HF1","HF2",
      "HF3","CL","HL1","HL2","HL3"],"Aliphatic")
]
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,label="Total",subplot=2)
for group in groups
    group_contrib = contrib(solute,mddf_dmf_acr.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=2)
end
plot!(
    xlim=(1,5),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF", subplot=2
)
savefig("./results/mddf_groups.png")

# 2D plot of group contributions
mers = collect(eachresidue(acr))
mer_contrib = zeros(length(mddf_dmf_acr.d),length(mers))
for (imer, mer) in pairs(mers)
    mer_contrib[:,imer] .= movavg(contrib(solute,mddf_dmf_acr.solute_atom,mer),n=10).x
end
idmin = findfirst( d -> d > 1.5, mddf_dmf_acr.d)
idmax = findfirst( d -> d > 3.5, mddf_dmf_acr.d)
labels = [ "$i" for i in 1:length(mers) ]
contourf(
    1:length(mers),
    mddf_dmf_acr.d[idmin:idmax],
    mer_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xlabel="mer",ylabel=L"r/\AA",
    xticks=(1:length(mers),labels),
)
savefig("./results/map2D_acr.png")








