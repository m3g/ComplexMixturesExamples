# ComplexMixtures.jl - Example

## POPC membrane solvated by water and ethanol

In this example ComplexMixtures.jl is used to study the interactions of a POPC membrane with a mixture of 20%(mol/mol) ethanol in water. At this concentration ethanol distabilizes the membrane. 

<center><img width=400px src="./system.png"></center>

System image: a POPC membrane (center) solvated by a mixture of water (purple) and ethanol (green). The system is composed by 59 POPC, 5000 water, and 1000 ethanol molecules.  

## Distribution functions and KB integrals 

![](./results/mddf_kb.png)

## Ethanol group contributions

![](./results/mddf_ethanol_groups.png)

## Interaction of POPC groups with water

![](./results/mddf_popc_water_groups.png)

## Interaction of POPC groups with ethanol

![](./results/mddf_popc_ethanol_groups.png)

## Ethanol interaction with the aliphatic chains

![](./results/map2D_aliphatic_chains.png)

## References

Membrane built with the VMD membrane plugin. 

Water and ethanol layers added with Packmol.

Density from: https://wissen.science-and-fun.de/chemistry/chemistry/density-tables/ethanol-water-mixtures/




