using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

# Plot defaults
function fig()

plot_font = "Computer Modern"
default(
    fontfamily=plot_font,
    linewidth=2.5, 
    framestyle=:box, 
    label=nothing, 
    grid=false,
    palette=:tab10
)
scalefontsizes() 
#scalefontsizes(1.3)

# Load system PDB file
system = readPDB("../simulation/equilibrated.pdb")

# Select the atoms corresponding to glycerol and water
dmf = select(system,"resname DMF")
acr = select(system,"resname FACR or resname ACR or resname LACR")

# Set ComplexMixtures selections
solute = Selection(acr,nmols=1)
solvent = Selection(dmf,natomspermol=12)

# Load previously computed data
mddf_dmf_acr = load("./mddf_dmf_acr.json")


l = @layout [ grid(2,1); grid(2,1); c{0.3h} ]
plot(layout=l)

# Plot the MDDF of DMF relative to PolyACR and its corresponding KB integral
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,ylabel="MDDF",subplot=1)
plot!(xlim=(0,20),subplot=1,xticks=:none)

y = movavg(mddf_dmf_acr.kb,n=10).x
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",subplot=2)
plot!(xlim=(0,20),subplot=2)

# Plot DMF group contributions to the MDDF. We use a vector of tuples here, 
# the first element of the tuple is a list of atom names, and the second element
# is the label of the group, to be used in the plot. 
groups = [ 
    (["C","O"],"CO"), # carbonyl
    (["HA"],"HA"),
    (["N"],"N"),
    (["CC","CT","HC1","HC2","HC3","HT1","HT2","HT3"],"Methyl groups")
]
x = mddf_dmf_acr.d
y = movavg(mddf_dmf_acr.mddf,n=10).x
plot!(x,y,label="Total",subplot=3,xticks=:none)
for group in groups
    # Retrieve the contributions of the atoms of this group
    group_contrib = contributions(solvent,mddf_dmf_acr.solvent_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=3)
end
plot!(
    xlim=(1,5),
    ylabel="MDDF",
    subplot=3
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
plot!(x,y,label="Total",subplot=4)
for group in groups
    group_contrib = contributions(solute,mddf_dmf_acr.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=4)
end
plot!(
    xlim=(1,5),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF", subplot=4
)

# 2D plot of group contributions
groups = [
    (["CF","HF1","HF2","HF3"],L"\textrm{CH_3}"),
    (["OE1","CD"],"CO"),
    (["NE2","HE22","HE21"],L"\textrm{NH_2}"),
    (["C","H2","H1","CA","HA"],L"\textrm{CHCH_2}"),
    (["CL","HL1","HL2","HL3"],L"\textrm{CH_3}")
]

# Here we split the polymer in residues, to extract the contribution of 
# each chemical group of each polymer mer independently
mers = collect(eachresidue(acr))
group_contribs = zeros(length(mddf_dmf_acr.d),0)
labels = String[]
for (imer, mer) in pairs(mers)
    # The mer.atoms array is just the original array of atoms,
    # and the atoms of this residue are retrieved by choosing the 
    # corresponding range (the eachresidue function does not copy)
    mer_atoms = mer.atoms[mer.range] 
    for igroup in 1:length(groups)
        group = groups[igroup]
        (imer != 1 && igroup == 1) && continue # only first residue has a first CH3
        (imer != 5 && igroup == 5) && continue # only last residue has a terminal CH3
        # And from the mer_atoms atoms, filter the ones corresponding to this group
        atoms = filter(at -> at.name in group[1], mer_atoms)
        # Retrive the contribution of these atoms to the MDDF
        contribs = movavg(contributions(solute,mddf_dmf_acr.solute_atom,atoms),n=10).x
        # Concatenate the results to build the 2D matrix
        group_contribs = hcat(group_contribs,contribs)
        # Push label to label list
        push!(labels,group[2])
    end
end

# Find the indices of the limits of the map we want
idmin = findfirst( d -> d > 1.5, mddf_dmf_acr.d)
idmax = findfirst( d -> d > 3.2, mddf_dmf_acr.d)

# Plot contour map
contourf!(
    1:length(labels),
    mddf_dmf_acr.d[idmin:idmax],
    group_contribs[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    xlabel="PolyACR group",ylabel=L"\textrm{Distance/\AA}",xrotation=60,
    xticks=(1:length(labels),labels),
    subplot=5,
)


for (y,lab) in [(8.2,"A)"),
                (7.1,"B)"),
                (5.4,"C)"),
                (4.4,"D)"),
                (2.4,"E)")]
    annotate!(
        -2.3, y,
        text(lab,plot_font,12),
        subplot=5,
    )
end

annotate!(1.1, 2.9, text("DMF\ngroups", plot_font, 10, :left), subplot=3)
annotate!(1.1, 2.9, text("PolyACR\ngroups", plot_font, 10, :left), subplot=4)

plot!(
    size=(550,800),
    leftmargin=10Plots.Measures.mm # adjust margin 
)

savefig("./polyacr.png")
savefig("./figure5.pdf")

end; fig()

