using ComplexMixtures
using PDBTools
using Plots
using LaTeXStrings
using EasyFit

function fig() # to simplify globals

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
scalefontsizes() 
#scalefontsizes(1.3)

# Load system PDB file
system = readPDB("../simulation/equilibrated.pdb")

# Select the atoms corresponding to glycerol and water
glyc = select(system,"resname GLLM")
water = select(system,"water")

# Set the solute structure for Glycerol, used to extract atomic contributions
solute = Selection(glyc,natomspermol=14)

# Load previously computed data
mddf_glyc = load("./mddf_glyc.json")
mddf_water_glyc = load("./mddf_water_glyc.json")

l = @layout [ 
    grid(2,1)
    grid(2,1)
    grid(2,1)
]
p = plot(layout=l)

# Plot the correlation functions
x = mddf_glyc.d # distances
y = movavg(mddf_glyc.mddf,n=10).x # the mddf (using movavg to smooth noise)
plot!(p,x,y,label="Glycerol-Glycerol",subplot=1)
x = mddf_water_glyc.d
y = movavg(mddf_water_glyc.mddf,n=10).x
plot!(p,x,y,label="Water-Glycerol",subplot=1)
plot!(p,
    ylabel="MDDF",
    xlim=(0,20),subplot=1,
    xticks=:none
)

# Plot the KB integrals
y = movavg(mddf_glyc.kb,n=10).x
plot!(p,x,y,subplot=2)
y = movavg(mddf_water_glyc.kb,n=10).x
plot!(p,x,y,subplot=2)
plot!(p,
    xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~cm^3~mol^{-1}}",
    xlim=(0,20),subplot=2
)

# Plot some group contributions to the MDDF. We select the atom names
# corresponding to each type of group of the glycerol molecule.  
hydroxyls = ["O1","O2","O3","HO1","HO2","HO3"]
aliphatic = ["C1","C2","C3","H11","H12","H2","H31","H32"]

# Retrieve the contributions of these groups to the MDDFs
hydroxyl_contrib = contrib(solute,mddf_glyc.solvent_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_glyc.solvent_atom,aliphatic)

# Plot group contributions
x = mddf_glyc.d
y = movavg(mddf_glyc.mddf,n=10).x
plot!(p,x,y,label="Total",subplot=3)
y = movavg(hydroxyl_contrib,n=10).x
plot!(p,x,y,label="Hydroxyl contributions",subplot=3)
y = movavg(aliphatic_contrib,n=10).x
plot!(p,x,y,label="Aliphatic contributions",subplot=3)
plot!(p,
    xlim=(1,6),
    ylim=(0,3.8),
    ylabel="MDDF",
    subplot=3,
    xticks=:none
)

x = mddf_water_glyc.d
y = movavg(mddf_water_glyc.mddf,n=10).x
plot!(p,x,y,subplot=4)
hydroxyl_contrib = contrib(solute,mddf_water_glyc.solute_atom,hydroxyls)
aliphatic_contrib = contrib(solute,mddf_water_glyc.solute_atom,aliphatic)

y = movavg(hydroxyl_contrib,n=10).x
plot!(p,x,y,subplot=4)
y = movavg(aliphatic_contrib,n=10).x
plot!(p,x,y,subplot=4)
plot!(p,
    xlim=(1,6),
    ylim=(0,3.8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
    subplot=4
)

# 2D maps plot of group contributions

# Glycerol groups
groups = [
    (["O1","HO1"],"OH"), # first hydroxyl
    (["C1","H11","H12"],L"\textrm{CH_2}"), # first CH2
    (["O2","HO2"],"OH"), # second hydroxyl
    (["C2","H2"],"CH"), # CH
    (["C3","H31","H32"],L"\textrm{CH_2}"), # second CH2
    (["O3","HO3"],"OH") # third hydroxyl
] 

# Contributions to Glycerol-Glycerol autocorrelation
group_contrib = zeros(length(mddf_glyc.d),length(groups))
for (igroup, group) in pairs(groups)
    group_contrib[:,igroup] .= contrib(solute,mddf_glyc.solute_atom,group[1])
end

idmin = findfirst( d -> d > 1.5, mddf_glyc.d)
idmax = findfirst( d -> d > 3.0, mddf_glyc.d)
labels = [ "OH", L"\textrm{CH_2}", "OH", "CH", L"\textrm{CH_2}", "OH" ] 

contourf!(p,
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xticks=:none,
    ylabel=L"\textrm{Distance/\AA}",
    subplot=5
)

# Water-glycerol interactions (Glycerol contributions)
group_contrib = zeros(length(mddf_glyc.d),length(groups))
for (igroup, group) in pairs(groups)
    group_contrib[:,igroup] .= contrib(solute,mddf_water_glyc.solute_atom,group[1])
end
contourf!(p,
    1:length(groups),
    mddf_glyc.d[idmin:idmax],
    group_contrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=5,
    xticks=(1:length(groups),labels),xrotation=60,
    ylabel=L"\textrm{Distance/\AA}",
    subplot=6
)
plot!(p,
    xlabel="Glycerol group",
    subplot=6
)


for (y,lab) in [(13.5,"A)"),
                (11.7,"B)"),
                (8.8,"C)"),
                (7.0,"D)"),
                (4.10,"E)"),
                (2.20,"F)")]
    annotate!(p,
        0, y,
        text(lab,plot_font,10),
        subplot=6,
    )
end

annotate!(p,1.1, 3.4, text("Glycerol-Glycerol", plot_font, 10, :left), subplot=3)
annotate!(p,1.1, 3.4, text("Water-Glycerol", plot_font, 10, :left), subplot=4)

annotate!(p, 4, 1.7, text("Glycerol-Glycerol", plot_font, 10, :left), subplot=5)
annotate!(p, 4, 1.7, text("Water-Glycerol", plot_font, 10, :left), subplot=6)

plot!(p,
    size=(550,800),
    leftmargin=10Plots.Measures.mm # adjust margin 
)

savefig(p,"./glyc_water.png")
savefig(p,"./figure3.pdf")

return p
end; 
fig()
