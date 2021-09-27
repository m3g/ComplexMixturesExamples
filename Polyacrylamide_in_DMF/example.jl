using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

# Plot defaults
plot_font = "Computer Modern"
default(fontfamily=plot_font,linewidth=2, framestyle=:box, label=nothing, grid=false)
scalefontsizes(1.3)

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
    (select(dmf,"name O"),"O"),
    (select(dmf,"name HA"),"HA"),
    (select(dmf,"name N"),"N"),
    (select(dmf,"name CC or 
                 name CT or 
                 name HC1 or 
                 name HC2 or
                 name HC3 or
                 name HT1 or
                 name HT2 or
                 name HT3"),"Methyl groups"),
]

x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot(x,y,label="Total")
for group in groups
    group_contrib = contrib(solvent,mddf_dmf_acr.solvent_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2])
end
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./results/mddf_dmf_groups.png")

# Plot ACR group contributions to the MDDF
groups = [ 
    (select(acr,"name OE1 or name CD"),"CO"),
    (select(acr,"name NE2 or
                 name HE22 or
                 name HE21"),L"\textrm{NH_2}"),
    (select(acr,"name C or
                 name H2 or
                 name H1 or
                 name CA or
                 name HA or
                 name CF or 
                 name HF1 or
                 name HF2 or
                 name HF3 or
                 name CL or
                 name HL1 or
                 name HL2 or
                 name HL3"),L"\textrm{Aliphatic}")
]

x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot(x,y,label="Total")
for group in groups
    group_contrib = contrib(solute,mddf_dmf_acr.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2])
end
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF"
)
savefig("./results/mddf_acr_groups.png")

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








