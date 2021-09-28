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

# Load system PDB file
system = readPDB("./simulation/equilibrated.pdb")

# Select the atoms corresponding to glycerol and water
glyc = select(system,"resname GLLM")
water = select(system,"water")

# Set the solute structure for Glycerol, used to extract atomic contributions
solute = Selection(glyc,natomspermol=14)

# Load previously computed data
mddf_glyc = load("./results/mddf_glyc.json")
mddf_water_glyc = load("./results/mddf_water_glyc.json")

# Plot the correlation functions and KB integrals
plot(layout=(2,1))
x = mddf_glyc.d
y = movavg(mddf_glyc.mddf,n=10).x
plot!(x,y,label="Glycerol")
x = mddf_water_glyc.d
y = movavg(mddf_water_glyc.mddf,n=10).x
plot!(x,y,label="Water")
plot!(
    ylabel="MDDF",
    xlim=(1.5,8),subplot=1
)

y = movavg(mddf_glyc.kb,n=10).x
plot!(x,y,subplot=2)
y = movavg(mddf_water_glyc.kb,n=10).x
plot!(x,y,subplot=2)
plot!(
    xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",
    xlim=(0,20),subplot=2
)
savefig("./results/mddf_kb.png")

# Plot some group contributions to the MDDF
hydroxyls = ["O1","O2","O3","HO1","HO2","HO3"]
aliphatic = ["C1","C2","C3","H11","H12","H2","H31","H32"]

hydroxyl_contrib = contrib(solute,mddf_glyc.solvent_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_glyc.solvent_atom,aliphatic)

plot(layout=(2,1))
x = mddf_glyc.d
y = movavg(mddf_glyc.mddf,n=10).x
plot!(x,y,label="Total",supblot=1)
y = movavg(hydroxyl_contrib,n=10).x
plot!(x,y,label="Hydroxyl contributions",subplot=1)
y = movavg(aliphatic_contrib,n=10).x
plot!(x,y,label="Aliphatic contributions",subplot=1)
plot!(
    xlim=(1,8),
    ylabel="MDDF",
    subplot=1
)

x = mddf_water_glyc.d
y = movavg(mddf_water_glyc.mddf,n=10).x
plot!(x,y,label="Total",subplot=2)
hydroxyl_contrib = contrib(solute,mddf_water_glyc.solute_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_water_glyc.solute_atom,aliphatic)

y = movavg(hydroxyl_contrib,n=10).x
plot!(x,y,label="Hydroxyl contributions",subplot=2)
y = movavg(aliphatic_contrib,n=10).x
plot!(x,y,label="Aliphatic contributions",subplot=2)
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
    subplot=2
)
savefig("./results/mddf_groups.png")

# 2D plot of group contributions
groups = [
    (["O1","HO1"],"OH"), # first hydroxyl
    (["C1","H11","H12"],L"\textrm{CH_2}"), # first CH2
    (["O2","HO2"],"OH"), # second hydroxyl
    (["C2","H2"],"CH"), # CH
    (["C3","H31","H32"],L"\textrm{CH_2}"), # second CH2
    (["O3","HO3"],"OH") # third hydroxyl
] 
group_contrib = zeros(length(mddf_glyc.d),length(groups))
for (igroup, group) in pairs(groups)
    group_contrib[:,igroup] .= contrib(solute,mddf_glyc.solute_atom,group[1])
end

idmin = findfirst( d -> d > 1.5, mddf_glyc.d)
idmax = findfirst( d -> d > 3.0, mddf_glyc.d)
labels = [ "OH", L"\textrm{CH_2}", "OH", "CH", L"\textrm{CH_2}", "OH" ] 

plot(layout=(2,1))
contourf!(
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xticks=(1:length(groups),labels),xrotation=60,
    ylabel=L"r/\AA",
    subplot=1
)

group_contrib = zeros(length(mddf_glyc.d),length(groups))
for (igroup, group) in pairs(groups)
    group_contrib[:,igroup] .= contrib(solute,mddf_water_glyc.solute_atom,group[1])
end
contourf!(
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xticks=(1:length(groups),labels),xrotation=60,
    ylabel=L"r/\AA",
    subplot=2
)
plot!(
    xlabel="Glycerol group",
    subplot=2
)

savefig("./results/map2D_glyc_glyc.png")
