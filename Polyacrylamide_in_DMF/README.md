# ComplexMixtures.jl - Example

## Polyacrylamide in DMF

In this example we illustrate how the solvation structure of a polymer can be studied with ComplexMixtures.jl. The system is a 5-mer segment of polyacrylamide (PAE - capped with methyl groups), solvated with dimethylformamide (DMF). The system is interesting because of the different functional groups and polarities involved in the interactions of DMF with PAE. A snapshot of the system is shown below. 

<center><img height=300px src="./system.png"></center>
The structures of DMF and of the polyacrylamide segment are:

<center>
<table><tr>
<td><img src=./simulation/dmf.png height=130px></td>
<td><img src=./simulation/polyacrylamide.png height=150px></td>
</tr>
<tr>
<td align=center>DMF</td>
<td align=center>Polyacrylamide</td>
</tr>
</table>
</center>

The system simulated consists of 2000 DMF molecules solvating a 5-mer chain of polyacrylamide. 

The step by step of this example is split into running the MDDF calculation, in the [compute_mddf.jl](./compute_mddf.jl) file, and extracting the information and plotting, in the [plots.jl](./plots.jl) file. 

The trajectory file, required to run the `compute_mddf.jl` script, is available [here - 275Mb](https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing), and should be saved
in the `test/trajectories` directory. The scripts

The `plots.jl` script can be executed from the results saved in this repository. 

To run these scripts directly, do:

```
cd ComplexMixturesExamples/Polyacrylamide_in_DMF
julia -e compute_mddf.jl
julia -e plots.jl
```

## Minimum-distance distribution function and KB integral

The distribution of DMF molecules around polyacrylamide is shown below. There is a peak at ~2.5Angs, indicating favorable non-specific interactions between the solvent molecules and the polymer. The peak is followed by a dip and diffuse peaks at higher distances. Thus, the DMF molecules are structured around the polymer, but essentially only in the first solvation shell.  

![](./results/mddf_kb.png)

The KB integral in a bicomponent mixture converges to the (negative of the) apparent molar volume of the solute. It is negative, indicating that the accumulation of DMF in the first solvation shell of the polymer is not enough to compensate the excluded volume of the solute. 

## Group contributions

The MDDF can be decomposed into the contributions of the DMF chemical groups, and on the polyacrylamide chemical groups. In the first panel below we show the contributions of the DMF chemical groups to the distribution function. 

The decomposition reveals that specific interactions peaking at distances slightly smaller than 2Angs exist between the polymer and the carbonyl group of DMF. Thus, there hydrogen bonds between the polymer and this group, which dominate the interactions between the solute and the solvent at short distances. The non-specific interactions peak at 2.5Angs and are composed of contributions of all DMF chemical groups, but particularly of the methyl groups. 

![](./results/mddf_groups.png)

The decomposition of the same MDDF in the contributions of the chemical groups of the polymer is clearly associated to the DMF contributions. The specific, hydrogen-bonding, interactions, are associated to the polymer amine groups. The amine groups also contribute to the non-specific interactions at greater distances, but these are a sum of the contributions of all polymer groups, polar or aliphatic. 

## Solvation along the polymer chain

We can decompose the MDDF into the contributions of each portion of the polymer chain. The map below displays the contributions of each chemical group of the polymer, now split into the mers of the polymer, to the MDDF.

![](./results/map2D_acr.png)

The terminal methyl groups interact strongly with DMF, and strong local density augmentations are visible in particular on the amine groups. These occur at less than 2.0Angs and are characteristic of hydrogen-bond interactions. Interestingly, the DMF molecules are excluded from the aliphatic and carbonyl groups of the polymer, relative to the other groups. 

Finally, it is noticeable that the central mer is more weakly solvated by DMF than the mers approaching the extremes of the polymer chain. This is likely a result of the partial folding of the polymer, that protects that central mers from the solvent in a fraction of the polymer configurations. 

## References

Molecules built with JSME: B. Bienfait and P. Ertl, JSME: a free molecule editor in JavaScript, Journal of Cheminformatics 5:24 (2013)
http://biomodel.uah.es/en/DIY/JSME/draw.en.htm

The system was built with [Packmol](http://m3g.iqm.unicamp.br/packmol).

The simulations were perfomed with [NAMD](https://www.ks.uiuc.edu/Research/namd/), with [CHARMM36](https://www.charmm.org) parameters. 
