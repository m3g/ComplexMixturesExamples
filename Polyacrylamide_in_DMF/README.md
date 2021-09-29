# ComplexMixtures.jl - Example

## Polyacrylamide in DMF

In this example we illustrate how the solvation structure of a polymer can be studied with ComplexMixtures.jl. The system is a 5-mer segment of polyacrylamide (PAE - capped with methyl groups), solvated with dimethylformamide (DMF). The system is interesting because of the different functional groups and polarities involved in the interactions of DMF with PAE. A snapshot of the system is shown below. 

<center><img height=300px src="./system.png"></center>

Trajectory file: [traj_Polyacry.dcd - 275Mb](https://drive.google.com/file/d/1ug43ncCLsBATaJrT9zlbaqK6AORVvhhx/view?usp=sharing)

The structures of DMF and of the polyacrylamide segment are:

<center>
<table><tr>
<td><img src=./simulation/dmf.png height=150px></td>
<td><img src=./simulation/polyacrylamide.png height=150px></td>
</tr>
<tr>
<td align=center>DMF</td>
<td align=center>Polyacrylamide</td>
</tr>
</table>
</center>

## Minimum-distance distribution function and KB integral

![](./results/mddf_kb.png)

## Group contributions

![](./results/mddf_groups.png)

## Solvation along the polymer chain

![](./results/map2D_acr.png)

## References

Molecules built with JSME:

B. Bienfait and P. Ertl, JSME: a free molecule editor in JavaScript, Journal of Cheminformatics 5:24 (2013)
http://biomodel.uah.es/en/DIY/JSME/draw.en.htm

System built with Packmol: 




