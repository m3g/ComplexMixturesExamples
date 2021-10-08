# ComplexMixtures.jl - Examples

This repository contains some detailed examples of the use of the [ComplexMixtures.jl](https://m3g.github.io/ComplexMixtures.jl/stable) package, for the study of the solvation structures of complex systems.

## Content

- [Protein solvated by water and glycerol](https://github.com/m3g/ComplexMixturesExamples/tree/main/Protein_in_Glycerol)

- [Crowded solution of glycerol in water](https://github.com/m3g/ComplexMixturesExamples/tree/main/Glycerol_in_Water)

- [Polyacrylamide in DMF](https://github.com/m3g/ComplexMixturesExamples/tree/main/Polyacrylamide_in_DMF)

- [POPC membrane in water/ethanol](https://github.com/m3g/ComplexMixturesExamples/tree/main/POPC_in_Water-Ethanol)

![image](./Protein_in_Glycerol/Density2D/density2D.png)

## New to Julia?

### Installing the package dependencies

The packages that are used in the script, with, for example, `using PDBTools`, can be installed in Julia using:

```julia-repl
julia> import Pkg

julia> Pkg.add("PDBTools")

```

or, if you simply type `]` your prompt will become the package manager prompt, `pkg>`, and then you do:

```
(@v1.6) pkg> add PDBTools, ComplexMixtures, LaTeXStrings

```

### Running the scripts

Julia compiles the code the first time it is executed. Thus, if you run one of the plotting scripts with, for example
```
julia plots.jl
```
it will take some time (a minute, perhaps), to produce the figure. If that is the only time you are running the script,
that is fine. However, if you want to modify the script and test alternative, that waiting is annoying. In that case, 
*start julia once*, and within `julia`, do:
```
julia> include("./plots.jl")
```

The first time that will take that annoying minute, but if you change the script and include the file again (press the up-arrow
to repeat the command), the new execution will be very fast. If you are into Julia and want to learn more sophisticated
workflows for a heavier development of code, see [this post](https://m3g.github.io/JuliaNotes.jl/stable/workflow/). 








