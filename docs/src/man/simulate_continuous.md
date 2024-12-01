```@setup tree_trait
using PhyloNetworks
using PhyloTraits
mkpath("../assets/figures")
truenet = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
```

# Continuous trait simulation

We show here how PhyloPlots can be used to simulate data under a Brownian motion
(BM) process. We use the same network as in the previous section.
The data used in the previous section was actually obtained using this
simulation (followed by some rounding).

## Trait simulation

We simulate three traits on the network: two independent traits that serve
as predictors, and one "dependent" trait that is affected by the first 2.
We start
by choosing the parameters of the BM (ancestral mean and variance), by creating
objects of class [`ParamsBM`](@ref)`<:ParamsProcess`.

```@example tree_trait
params_trait1 = ParamsBM( 2, 0.5) # BM with mean  2 and variance 0.5
params_trait2 = ParamsBM(-2, 1)   # BM with mean -2 and variance 1.0
nothing # hide
```

We then simulate the traits to serve as predictors according to these parameters, using
function [`rand`](@ref).
For reproducibility, we use a stable random number generator (RNG),
but the default random generator is better (and more efficient).
To use the default RNG, simply remove `rng` in the code below.

```@example tree_trait
using StableRNGs; rng = StableRNG(18480224); # stable RNG for reproducibility
sim1 = rand(rng, truenet, params_trait1)  # simulate a BM on truenet
sim2 = rand(rng, truenet, params_trait2)
nothing # hide
```

This creates objects of class [`TraitSimulation`](@ref), from which we can
extract the data at the tips, thanks to the method
[`getindex(::TraitSimulation, ::Symbol)`](@ref).

```@example tree_trait
trait1 = sim1[:tips] # trait 1 at the tips (data)
trait2 = sim2[:tips]
nothing # hide
```

This extractor creates an `Array` with one column, and as many lines as the
number of tips there are in the phylogeny.  It is sorted in the same order as
the tips of the phylogeny used to simulate it.  
If needed, we could also extract the simulated values at the internal nodes
in the network:

```@example tree_trait
sim1[:internalnodes]
nothing # hide
```
These values are those used in the previous section, when we compared the
ancestral state reconstruction (or "predictions") of trait 1 (stored in
`ancTrait1`) to the true simulated value.

Finally, we generate the last trait correlated with trait 1
(but not trait 2), with phylogenetic noise.
```@example tree_trait
noise = rand(StableRNG(18700904), truenet, ParamsBM(0, 0.1)) # phylogenetic residuals
trait3 = 10 .+ 2 * trait1 .+ noise[:tips] # trait to study. independent of trait2
nothing # hide
```

For a simulation study, we would assume that we measured the three traits above,
and then we would like to apply analysis methods to these data, e.g. to
see how they infer the impact of traits 1 and 2 on trait 3.
To do that, we may need to create a data frame containing all 3 traits,
as shown below:

```@repl tree_trait
using DataFrames
dat = DataFrame(trait1 = trait1, trait2 = trait2, trait3 = trait3,
                tipnames = tiplabels(sim1))
```

These data were rounded to 3 digits and used in the previous section.
