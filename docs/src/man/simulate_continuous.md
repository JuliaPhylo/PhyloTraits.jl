```@setup sim_BM
using PhyloNetworks
using PhyloTraits
using PhyloPlots, RCall
mkpath("../assets/figures")
figname(x) = joinpath("..", "assets", "figures", x)
truenet = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
```

# Continuous trait simulation

We show here how PhyloPlots can be used to simulate data under a Brownian motion
(BM) process. We use the same network as in the previous section.
The data used in the previous section was actually obtained using this
simulation (followed by some rounding).

## Brownian motion

We simulate three traits on the network: two independent traits that serve
as predictors, and one "dependent" trait that is affected by the first 2.
We start
by choosing the parameters of the BM (ancestral mean and variance), by creating
objects of class [`ParamsBM`](@ref)`<:ParamsProcess`.

```@example sim_BM
params_trait1 = ParamsBM( 2, 0.5) # BM with mean  2 and variance 0.5
params_trait2 = ParamsBM(-2, 1)   # BM with mean -2 and variance 1.0
nothing # hide
```

We then simulate the traits to serve as predictors according to these parameters, using
function [`rand`](@ref).
For reproducibility, we use a stable random number generator (RNG),
but the default random generator is better (and more efficient).
To use the default RNG, simply remove `rng` in the code below.

```@example sim_BM
using StableRNGs; rng = StableRNG(18480224); # stable RNG for reproducibility
sim1 = rand(rng, truenet, params_trait1)  # simulate a BM on truenet
sim2 = rand(rng, truenet, params_trait2)
nothing # hide
```

This creates objects of class [`TraitSimulation`](@ref), from which we can
extract the data at the tips, thanks to the method
[`getindex(::TraitSimulation, ::Symbol)`](@ref).

```@repl sim_BM
trait1 = sim1[:tips] # trait 1 at the tips (data)
trait2 = sim2[:tips]
```

This extractor creates an `Array` with one column, and as many lines as the
number of tips there are in the phylogeny.  It is sorted in the same order as
the tips of the phylogeny used to simulate it.  
If needed, we could also extract the simulated values at the internal nodes
in the network:

```@example sim_BM
sim1[:internalnodes]
```
These values (rounded to 3 digits) are those used in section
[Ancestral state reconstruction](@ref),
when we compared the ancestral state reconstruction (or "predictions") of
trait 1 (stored in `ancTrait1`) to the true simulated value.

Finally, we generate the last trait correlated with trait 1
(but not trait 2), with phylogenetic noise.
```@example sim_BM
noise = rand(StableRNG(18700904), truenet, ParamsBM(0, 0.1)) # phylogenetic residuals
trait3 = 10 .+ 2 * trait1 .+ noise[:tips] # trait to study. independent of trait2
nothing # hide
```

For a simulation study, we would assume that we measured the three traits above,
and then we would like to apply analysis methods to these data, e.g. to
see how they infer the impact of traits 1 and 2 on trait 3.
To do that, we may need to create a data frame containing all 3 traits,
as shown below:

```@repl sim_BM
using DataFrames
dat = DataFrame(trait1 = trait1, trait2 = trait2, trait3 = trait3,
                tipnames = tiplabels(sim1))
```

These data were rounded to 3 digits and used in the previous section on
[Phylogenetic regression](@ref).

## shifted Brownian motion

In a shifted BM, the trait evolves as a BM on the network most of
the time, but undergo sudden *shifts* on some of the branches.
The positions and values of the shifts can be stored in a [`ShiftNet`](@ref)
object. The position of the shifts can be given using vector of edges.
To see this, let's first plot the network with its associated edges and node
numbers.
```@example sim_BM
R"svg"(figname("truenet_with_numbers.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, useedgelength=true, showedgenumber=true);
R"dev.off()" # hide
nothing # hide
```
![truenet_with_numbers](../assets/figures/truenet_with_numbers.svg)

For identifiability reasons, shifts are only allowed on tree (not hybrid) edges.
Here, a shift on either hybrid edge 9 or 7 would have the same
effect as a shift on edge 6: by shifting traits of species A and B.

Let's say that we want to add a shift with value 5.0 on the branch directly
following the hybridization event, in order to model transgressive evolution.
We can see on the
plot above that this branch is number 6, so we define the following object:

```@repl sim_BM
edge_afterhyb = truenet.edge[6] # the 6th edge has number 6
shift = ShiftNet(truenet.edge[6], 5.0,  truenet)
```
Note that the edge numbers and values of a `ShiftNet` object can be retrieved
thanks to functions [`getshiftedgenumber`](@ref) and [`getshiftvalue`](@ref).
The constructor can take a single edge and associated value, like here,
or two vectors of edges and matching values.

Because we often need to put shifts only on edges right after hybrids,
there is a special function [`shiftathybrids`](@ref) to do that, so that
we do not have to find out their edges number. Here, the `shift` object
could hence have been defined as:
```@example sim_BM
shift = shiftathybrids(5.0,  truenet)
```

The parameters for the simulation are then defined as above, just adding
our `ShiftNet` object as a parameter.
Again, we use a custom stable random number generator to make our example
reproducible, but we recommend *against* this random number generator forf
simulation studies (just remove the argument `rng`).

```@example sim_BM
params_sh = ParamsBM(2, 0.5, shift) # BM with mean 2, variance 0.5, and shifts.
nothing # hide
```

The traits are simulated using the same function [`rand`](@ref), and
extracted at the tips as before.
```@example sim_BM
rng = StableRNG(18700904)
sim_sh = rand(rng, truenet, params_sh) # simulate a shifted BM on truenet
trait_sh = sim_sh[:tips]          # trait at the tips (data)
```
```@example sim_BM
tiplabels(sim_sh)
```
