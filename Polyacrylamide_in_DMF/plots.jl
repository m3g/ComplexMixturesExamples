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
dmf = select(system,"resname DMF")
acr = select(system,"resname FACR or resname ACR or resname LACR")

# Set ComplexMixtures selections
solute = Selection(acr,nmols=1)
solvent = Selection(dmf,natomspermol=12)

# Load previously computed data
mddf_dmf_acr = load("./results/mddf_dmf_acr.json")

# Plot the MDDF of DMF relative to PolyACR and its corresponding KB integral
plot(layout=(2,1))
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,ylabel="MDDF",subplot=1)
plot!(xlim=(0,20),subplot=1)

y = movavg(mddf_dmf_acr.kb,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",subplot=2)
plot!(xlim=(0,20),subplot=2)
savefig("./results/mddf_kb.png")

# Plot DMF group contributions to the MDDF
groups = [ 
    (["C","O"],"CO"), # carbonyl
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
    (["CF","HF1","HF2","HF3",
      "CL","HL1","HL2","HL3"],L"\textrm{CH_3}"), # terminal methyles
    (["OE1","CD"],"CO"), # carbonyl
    (["NE2","HE22","HE21"],L"\textrm{NH_2}"), # amine
    (["C","H2","H1","CA","HA"],L"\textrm{CHCH_2}"), # backbone
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
groups = [
    (["CF","HF1","HF2","HF3"],L"\textrm{CH_3}"),
    (["OE1","CD"],"CO"),
    (["NE2","HE22","HE21"],L"\textrm{NH_2}"),
    (["C","H2","H1","CA","HA"],L"\textrm{CHCH_2}"),
    (["CL","HL1","HL2","HL3"],L"\textrm{CH_3}")
]

mers = collect(eachresidue(acr))
group_contribs = zeros(length(mddf_dmf_acr.d),0)
labels = String[]
for (imer, mer) in pairs(mers)
    mer_atoms = mer.atoms[mer.range]
    for igroup in 1:length(groups)
        group = groups[igroup]
        (imer != 1 && igroup == 1) && continue
        (imer != 5 && igroup == 5) && continue
        atoms = filter(at -> at.name in group[1], mer_atoms)
        contribs = movavg(contrib(solute,mddf_dmf_acr.solute_atom,atoms),n=10).x
        group_contribs = hcat(group_contribs,contribs)
        push!(labels,group[2])
    end
end

idmin = findfirst( d -> d > 1.5, mddf_dmf_acr.d)
idmax = findfirst( d -> d > 3.2, mddf_dmf_acr.d)
contourf(
    1:length(labels),
    mddf_dmf_acr.d[idmin:idmax],
    group_contribs[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    xlabel="Group",ylabel=L"r/\AA",xrotation=60,
    xticks=(1:length(labels),labels),
    margin=5Plots.Measures.mm
)
savefig("./results/map2D_acr.png")






















