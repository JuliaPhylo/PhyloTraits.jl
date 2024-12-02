```@setup sim_discrete
using PhyloNetworks
using PhyloTraits
using DataFrames
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
```

# Discrete trait simulation

[`rand`](@ref) can simulate discrete traits along a known network.

## binary trait example

For example, we can define a binary trait model with states
"carnivory" (state 1) and "non-carnivory" (state 2), then ask for
a trait to be simulated along our network. We can ask for
3 independent simulations, giving us 3 traits then, arranged in 3 rows.

```@repl sim_discrete
m1 = BinaryTraitSubstitutionModel(1.0,2.0, ["carnivory", "non-carnivory"])
```

We also need a phylogeny. We re-use the network from the section on
[Discrete trait evolution](@ref).
We then simulate traits with a stable random number generator (RNG) for
reproducibility, but recommend using the default RNG by simply removing
the argument `rng` below.

```@repl sim_discrete
net = readnewick("(O:4,(A:3,((B:0.4)#H1:1.6::0.92,((#H1:0::0.08,C:0.4):0.6,(D:.2,E:.2):0.8):1):1):1);");
using StableRNGs; rng = StableRNG(123); # for reproducibility of this example
traitmatrix, nodecolumn = rand(rng, m1, net; ntraits=3);
traitmatrix
```

In this trait matrix, each column corresponds to a node,
each row is a trait, and each entry gives the state of that trait for that node,
as an index. To get the state labels:

```@repl sim_discrete
labs = getlabels(m1)
labs[traitmatrix]
```

The `nodecolumn` vector says which node corresponds to which column
in the trait matrix, and we can compare to the node numbers in the network.
For example, the first column corresponds to node `-2`, which is the root.
(The root is always in the first column: that's where the simulation starts.)

```@repl sim_discrete
nodecolumn
getroot(net)
```

As another example, we can extract the data for species "A".
We first get the column number for "A" (column 12),
then get A's data in that column.

```@repl sim_discrete
findfirst(isequal("A"), nodecolumn) # 12
nodecolumn[12]
traitmatrix[:,12] # data for species A, as category indices
labs[traitmatrix[:,12]] # same data, as category names
```

## example of DNA simulation

Below we simulate 4 sites of a DNA alignment, independent from an HKY model
with transition/transversion ratio κ=0.5 and stationary base frequencies of
0.2 for A and T and 0.3 for C and G:

and 
```@repl sim_discrete
m2 = HKY85([.5], [0.20, 0.30, 0.30, 0.20])
rng = StableRNG(36154); # again, for reproducibility
traitmatrix, nodecolumn = rand(rng, m2, net; ntraits=4);
traitmatrix
labs = getlabels(m2)
```

To get the data at the tips only, and in a specific order, we can do this.
```@repl sim_discrete
taxa = tiplabels(net)
taxoncol = indexin(taxa, nodecolumn) # column indices to get taxa in order
labs[traitmatrix[:,taxoncol]] # trait at the tips only, ordered as in 'taxa'
```

This type of DNA data is from the package
[BioSymbols](https://biojulia.dev/BioSymbols.jl/stable/),
for interoperability with packages from [BioJulia](https://biojulia.dev).

```@repl sim_discrete
using BioSequences
d = Dict(taxa[i] => LongDNA{4}(labs[traitmatrix[:,taxoncol[i]]])
         for i in eachindex(taxa)) 
```
