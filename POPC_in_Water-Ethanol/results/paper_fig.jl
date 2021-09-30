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
scalefontsizes()
#scalefontsizes(1.3)

# Load pdb file of the system
system = readPDB("../simulation/equilibrated.pdb")

# Select the atoms corresponding to glycerol and water
popc = select(system,"resname POPC")
water = select(system,"water")
ethanol = select(system,"resname ETOH")

# Set the complete membrane as a single solute
solute = Selection(popc,nmols=1)

# Load previously computed data
mddf_water_POPC = load("./mddf_water_POPC.json")
mddf_ethanol_POPC = load("./mddf_ethanol_POPC.json")

# Set a solvent structure, uses to retrieve the atomic Contributions
# of ethanol to the MDDF
solvent = Selection(ethanol,natomspermol=9)

l = @layout [ grid(2,1)  b{0.3w} ; 
              c          d{0.5w} ;
              grid(2,1) grid(2,1) ]
plot(layout=l)

# Water-POPC MDDF and KB integral. Here we use `movavg` from `EasyFit` to soften
# the noise
x = mddf_water_POPC.d # distances
y = movavg(mddf_water_POPC.mddf,n=10).x
plot!(x,y,label="Water",
    #xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",subplot=1
)

y = movavg(mddf_water_POPC.kb,n=10).x/1000
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~L~mol^{-1}}",subplot=2)
plot!(xlim=(0,10),subplot=2)

# Ethanol-POPC MDDF and KB integral
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot!(x,y,label="Ethanol")
plot!(xlim=(0,10),xticks=:none,subplot=1)

y = movavg(mddf_ethanol_POPC.kb,n=10).x/1000
plot!(x,y,xlabel=L"\textrm{Distance / \AA}",ylabel=L"\textrm{KB~/~L~mol^{-1}}",subplot=2)
plot!(xlim=(0,10),subplot=2)

# Contributions of the ethanol groups
groups = [
    (select(ethanol,"name O1 or name HO1"),"Hydroxyl"),
    (select(ethanol,"not name O1 and not name HO1"),"Aliphatic chain"),
]
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot!(x,y,label="Total",subplot=3)
for group in groups
    group_contrib = contrib(solvent,mddf_ethanol_POPC.solvent_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=3)
end
plot!(
    xlim=(1,8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF", subplot=3
)

# Contributions of POPC groups to ethanol-POPC distribution. Bellow is the list
# of the atom types of each POPC group. We need to list those to select specifically
# the atoms of each chemical group of interest. 
groups = [
  (["N","C12","H12A","H12B","C13","H13A","H13B","H13C","C14",
    "H14A","H14B","H14C","C15","H15A","H15B","H15C","C11","H11A","H11B"],"Choline"),
  ([ "P","O13","O14","O12" ],"Phosphate"),
  (["O11","C1","HA","HB","C2","HS","O21","C3","HX","HY","O31"],"Glycerol"),
  (["O22","C21","H2R","H2S","C22","C23","H3R","H3S","C24","H4R","H4S",
    "C25","H5R","H5S","C26","H6R","H6S","C27","H7R","H7S","C28","H8R","H8S",
    "C29","H91","C210","H101","C211","H11R","H11S","C212","H12R","H12S",
    "C213","H13R","H13S","C214","H14R","H14S","C215","H15R","H15S",
    "C216","H16R","H16S","C217","H17R","H17S","C218","H18R","H18S","H18T"],"Oleoyl"),
  (["C31","O32","C32","H2X","H2Y","C33","H3X","H3Y","C34","H4X","H4Y",
    "C35","H5X","H5Y","C36","H6X","H6Y","C37","H7X","H7Y","C38","H8X",
    "H8Y","C39","H9X","H9Y","C310","H10X","H10Y","C311","H11X","H11Y",
    "C312","H12X","H12Y","C313","H13X","H13Y","C314","H14X","H14Y","C315",
    "H15X","H15Y","C316","H16X","H16Y","H16Z"],"Palmitoyl")
]
x = mddf_ethanol_POPC.d
y = movavg(mddf_ethanol_POPC.mddf,n=10).x
plot!(x,y,label="Total",subplot=4)
for group in groups
    # Retrieve the contributions of the atoms of this group
    group_contrib = contrib(solute,mddf_ethanol_POPC.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=4)
end
plot!(
    xlim=(1.3,5),
    ylim=(0,1.8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
    subplot=4
)

# the same for POPC-water
x = mddf_water_POPC.d
y = movavg(mddf_water_POPC.mddf,n=10).x
plot!(x,y,label="Total",subplot=5)
for group in groups
    group_contrib = contrib(solute,mddf_water_POPC.solute_atom,group[1])
    y = movavg(group_contrib,n=10).x
    plot!(x,y,label=group[2],subplot=5)
end
plot!(
    xlim=(1.3,5),
    ylim=(0,1.8),
    xlabel=L"\textrm{Distance / \AA}",
    ylabel="MDDF",
    subplot=5
)

# Map of interactions of ethanol with Oleoyl groups. Now we split the
# Oleoyl chain into its components, ordered along the chain. 
oleoyl_groups = [
  (["O22","C21"],"CO"),
  (["H2R","H2S","C22"],L"\textrm{CH_2}"),
  (["C23","H3R","H3S"],L"\textrm{CH_2}"),
  (["C24","H4R","H4S"],L"\textrm{CH_2}"),
  (["C25","H5R","H5S"],L"\textrm{CH_2}"),
  (["C26","H6R","H6S"],L"\textrm{CH_2}"),
  (["C27","H7R","H7S"],L"\textrm{CH_2}"),
  (["C28","H8R","H8S"],L"\textrm{CH_2}"),
  (["C29","H91"],L"\textrm{CH}"),
  (["C210","H101"],L"\textrm{CH}"),
  (["C211","H11R","H11S"],L"\textrm{CH_2}"),
  (["C212","H12R","H12S"],L"\textrm{CH_2}"),
  (["C213","H13R","H13S"],L"\textrm{CH_2}"),
  (["C214","H14R","H14S"],L"\textrm{CH_2}"),
  (["C215","H15R","H15S"],L"\textrm{CH_2}"),
  (["C216","H16R","H16S"],L"\textrm{CH_2}"),
  (["C217","H17R","H17S"],L"\textrm{CH_2}"),
  (["C218","H18R","H18S","H18T"],L"\textrm{CH_3}"),
]

gcontrib = zeros(length(mddf_ethanol_POPC.d),length(oleoyl_groups))
for (ig, group) in pairs(oleoyl_groups)
    gcontrib[:,ig] .= movavg(contrib(solute,mddf_ethanol_POPC.solute_atom,group[1]),n=10).x
end
labels = [ oleoyl_groups[i][2] for i in 1:length(oleoyl_groups) ]
idmin = findfirst( d -> d > 1.5, mddf_ethanol_POPC.d)
idmax = findfirst( d -> d > 3.0, mddf_ethanol_POPC.d)

contourf!(
    1:length(oleoyl_groups),
    mddf_ethanol_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    ylabel=L"\textrm{Distance/\AA}",xrotation=60,
    xticks=(1:length(oleoyl_groups),labels),
    subplot=6
)

# And with Palmitoyl groups
palmitoyl_groups = [
  (["C31","O32"],"CO"),
  (["C32","H2X","H2Y"],L"\textrm{CH_2}"),
  (["C33","H3X","H3Y"],L"\textrm{CH_2}"),
  (["C34","H4X","H4Y"],L"\textrm{CH_2}"),
  (["C35","H5X","H5Y"],L"\textrm{CH_2}"),
  (["C36","H6X","H6Y"],L"\textrm{CH_2}"),
  (["C37","H7X","H7Y"],L"\textrm{CH_2}"),
  (["C38","H8X","H8Y"],L"\textrm{CH_2}"),
  (["C39","H9X","H9Y"],L"\textrm{CH_2}"),
  (["C310","H10X","H10Y"],L"\textrm{CH_2}"),
  (["C311","H11X","H11Y"],L"\textrm{CH_2}"),
  (["C312","H12X","H12Y"],L"\textrm{CH_2}"),
  (["C313","H13X","H13Y"],L"\textrm{CH_2}"),
  (["C314","H14X","H14Y"],L"\textrm{CH_2}"),
  (["C315","H15X","H15Y"],L"\textrm{CH_2}"),
  (["C316","H16X","H16Y","H16Z"],L"\textrm{CH_3}")
]
gcontrib = zeros(length(mddf_ethanol_POPC.d),length(palmitoyl_groups))
for (ig, group) in pairs(palmitoyl_groups)
    gcontrib[:,ig] .= movavg(contrib(solute,mddf_ethanol_POPC.solute_atom,group[1]),n=10).x
end
labels = [ palmitoyl_groups[i][2] for i in 1:length(palmitoyl_groups) ]
idmin = findfirst( d -> d > 1.5, mddf_ethanol_POPC.d)
idmax = findfirst( d -> d > 3.0, mddf_ethanol_POPC.d)

# This creates a contour plot, with some options to make the plot look nicer.
contourf!(
    1:length(palmitoyl_groups),
    mddf_ethanol_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    xlabel="Group",ylabel=L"\textrm{Distance/\AA}",xrotation=60,
    xticks=(1:length(palmitoyl_groups),labels),subplot=7
)

annotate!( 10, 2.7, text("Ethanol-Oleoyl", :left, 12, plot_font), subplot=6)
annotate!(  8, 2.7, text("Ethanol-Palmitoyl", :left, 12, plot_font), subplot=7)

# And now the same for water-POPC

gcontrib = zeros(length(mddf_water_POPC.d),length(oleoyl_groups))
for (ig, group) in pairs(oleoyl_groups)
    gcontrib[:,ig] .= movavg(contrib(solute,mddf_water_POPC.solute_atom,group[1]),n=10).x
end
labels = [ oleoyl_groups[i][2] for i in 1:length(oleoyl_groups) ]
idmin = findfirst( d -> d > 1.5, mddf_water_POPC.d)
idmax = findfirst( d -> d > 3.0, mddf_water_POPC.d)

contourf!(
    1:length(oleoyl_groups),
    mddf_water_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    ylabel=L"\textrm{Distance/\AA}",xrotation=60,
    xticks=(1:length(oleoyl_groups),labels),subplot=8
)

gcontrib = zeros(length(mddf_water_POPC.d),length(palmitoyl_groups))
for (ig, group) in pairs(palmitoyl_groups)
    gcontrib[:,ig] .= movavg(contrib(solute,mddf_water_POPC.solute_atom,group[1]),n=10).x
end
labels = [ palmitoyl_groups[i][2] for i in 1:length(palmitoyl_groups) ]
idmin = findfirst( d -> d > 1.5, mddf_water_POPC.d)
idmax = findfirst( d -> d > 3.0, mddf_water_POPC.d)

contourf!(
    1:length(palmitoyl_groups),
    mddf_water_POPC.d[idmin:idmax],
    gcontrib[idmin:idmax,:],
    color=cgrad(:tempo),linewidth=1,linecolor=:black,
    colorbar=:none,levels=10,
    xlabel="Group",ylabel=L"\textrm{Distance/\AA}",xrotation=60,
    xticks=(1:length(palmitoyl_groups),labels),subplot=9
)

annotate!( 10, 2.7, text("Water-Oleoyl", :left, 12, plot_font), subplot=8)
annotate!(  8, 2.7, text("Water-Palmitoyl", :left, 12, plot_font), subplot=9)

for (c,lab) in [((-21.8, 16.4),"A)"),
                ((-21.8, 14.1),"B)"),
                (( -2.4, 15.2),"C)"),
                ((-21.8,  9.4),"D)"),
                (( -2.4,  9.4),"E)"),
                ((-21.8,  6.0),"F)"),
                (( -2.4,  6.0),"G)")]
    annotate!(
        c[1],c[2],
        text(lab,plot_font,12),
        subplot=9,
    )
end

annotate!(
    5,2,
    text("Ethanol contributions",plot_font,10),
    subplot=3
)
annotate!(
    3.8,0.3,
    text("PolyACR contributions\nto Ethanol MDDF",plot_font,10),
    subplot=4
)
annotate!(
    3.8,0.4,
    text("PolyACR contributions\nto Water MDDF",plot_font,10),
    subplot=5
)

plot!(
    leftmargin=10Plots.Measures.mm,
    size=(1000,800)
)

savefig("./popc.png")
