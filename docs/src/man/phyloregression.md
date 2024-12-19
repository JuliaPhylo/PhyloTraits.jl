```@setup tree_trait
using PhyloNetworks
using PhyloTraits
mkpath("../assets/figures")
```

# Continuous trait analysis

After inferring a phylogeny, we can take
these phylogenetic relationships into account when studying the distribution of
quantitative traits measured for extant species.
This is the goal of phylogenetic comparative methods (PCM).
With PhyloTraits, we can do so for a phylogeny that is either a tree,
or a network with reticulations.
More details can be found on the developments below in Bastide et al. 2018 [^B18]

We assume a fixed network (which may be a tree), correctly rooted, with branch
lengths proportional to calendar time.
Below, we use a network that is time-consistent (all paths from the root to any
given node have the same length) and ultrametric (all the tips are contemporary).

```@example tree_trait
truenet = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");
```
We can plot the network thanks to package
[PhyloPlots](https://github.com/juliaphylo/PhyloPlots.jl), which uses
[`RCall`](https://juliainterop.github.io/RCall.jl/stable/gettingstarted/).
The `name` function is only instrumental here, to ensure that the figure is
saved in the correct directory when the documentation is built.
We only show the commands to actually save the plot in this first example for
the interested reader, but we will hide those in the rest of the chapter, for
the sake of clarity.
```@example tree_trait
using PhyloPlots, RCall
name(x) = joinpath("..", "assets", "figures", x)
R"svg"(name("truenet.svg"), width=8, height=4)
R"par"(mar=[0,0,0,0])
plot(truenet, useedgelength=true, showgamma=true);
R"dev.off()"
nothing # hide
```
![truenet](../assets/figures/truenet.svg)

## Model and variance matrix

Assuming that the network is known and that the continuous traits evolve under a
Brownian Motion (BM) over time, it is possible to compute the expected variance
covariance matrix between tip measurements. This can be done using function
[`PhyloNetworks.vcv`](@extref), whose syntax is inspired from the well known
corresponding [`ape`](https://CRAN.R-project.org/package=ape) function.
```@repl tree_trait
C = vcv(truenet)
```
The matrix is returned as a `DataFrame`, with columns named by the
tips of the network to allow for easy identification.
Each row also corresponds to a tip in the network, and rows are
ordered in the same way as columns.

The computation of this matrix is based on the more general function
[` PhyloNetworks.sharedpathmatrix`](@extref).
It is at the core of all the Phylogenetic Comparative Methods described below.

## Phylogenetic regression

Assume that we measured three continuous traits in the data frame below.
We want to study the impact of traits 1 and 2 on trait 3.
To do that, we can perform a phylogenetic regression.

To make sure that data are mapped to the correct tip in the phylogeny,
even though the tips might be ordered differently in the data frame compared
to the phylogeny, our data needs to have a column with the names of the tips in
the network.
If this column is labeled `tipnames`, fitting the data will not require an
extra option.
```@example tree_trait
using DataFrames
dat = DataFrame(
  trait1 = [ 2.668,  3.696,  4.541, 4.846,  2.268, -0.331],
  trait2 = [-3.008, -4.146, -2.338, 0.655, -3.339, -4.566],
  trait3 = [15.424, 17.333, 18.115, 18.81, 13.337, 10.012],
  tipnames = ["D", "C", "A", "B", "E", "O"]
)
nothing # hide
```

Phylogenetic regression / ANOVA is based on the
[GLM](https://github.com/JuliaStats/GLM.jl) package, with the network as an
extra argument, using function [`phylolm`](@ref).

```@example tree_trait
using StatsModels # for statistical model formulas
fitTrait3 = phylolm(@formula(trait3 ~ trait1 + trait2), dat, truenet)
```
The REML criterion is used by default, for estimating the variance
parameter(s). ML could be used instead with option `reml=false`.  
From this, we can see that the intercept, the coefficient for trait 1
and the variance of the noise are correctly estimated
(given that there are only 6 taxa).
In addition, the Student T test for the coefficient
associated with trait 2 has a high p-value, which means that this coefficient
is not significantly different from 0. This is consistent with the
way we simulated trait 3.

The function returns an object of type [`PhyloNetworkLinearModel`](@ref)`<:GLM.LinPredModel`.
It is a subtype of the GLM type `LinPredModel`, which means that all base
functions from Julia [StatsBase](https://github.com/JuliaStats/StatsBase.jl) can
be applied to it. See the documentation for this type for a list of all
functions that can be used. Some functions allow the user to retrieve directly
the estimated parameters of the BM, and are specific to this object.
```@repl tree_trait
sigma2_phylo(fitTrait3) # estimated variance of the BM
mu_phylo(fitTrait3) # estimated root value of the BM
```

## Ancestral state reconstruction

### From known parameters

If we assume that we know the exact model of evolution that generated the traits,
we can do ancestral trait reconstruction. Here, we simulated trait 1 ourselves
(see next section), so we can use the true process with the true parameters.
In other words, we can reconstruct the state at the internal nodes,
given the values at the tips, the known value at the root (2)
and the known BM variance (0.5).
```@example tree_trait
ancTrait1 = ancestralreconstruction(truenet, dat.trait1, ParamsBM(2, 0.5))
nothing # hide
```
Function [`ancestralreconstruction`](@ref) creates an object with type
[`ReconstructedStates`](@ref). Several extractors can be applied to it:
```@repl tree_trait
using StatsAPI, StatsBase # for predict and stderror()
predict(ancTrait1) # predictions
stderror(ancTrait1) # associated standard errors
predict(ancTrait1, interval=:prediction, level=0.90) # prediction interval (at level 90%)
```
We can plot the ancestral states or prediction intervals on the tree, using the
`nodelabel` argument of the `plot` function.
```@example tree_trait
R"svg"(name("ancestral_expe.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel=predict(ancTrait1, text=true), tipoffset=0.1);
R"dev.off()" # hide
nothing # hide
```
![ancestral_expe](../assets/figures/ancestral_expe.svg)

```@example tree_trait
ancInt = predict(ancTrait1, interval=:prediction, text=true) # format the prediction intervals for the plot
R"svg"(name("ancestral_predint.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel=ancInt[!,[:nodenumber,:interval]], tipoffset=0.1);
R"dev.off()" # hide
nothing # hide
```
![ancestral_predint](../assets/figures/ancestral_predint.svg)

The `predict` function has an optional argument to state
the `level` of the prediction interval. If not given, the default value is
0.95.

It is also possible to plot both the reconstructed state and the predicted value
on the same plot, using the optional keyword argument `combine`.
As shown below, we could also use the `RCall` method from the
[`plot`](https://juliaphylo.github.io/PhyloPlots.jl/stable/lib/public/) function.
```@example tree_trait
ancInt = predict(ancTrait1, interval=:prediction, text=true, combine=true) 
plot(truenet, nodelabel = ancInt[!,[:nodenumber,:interval]], tipoffset=0.1);
nothing # hide
```
These plots tend to be quite busy, even for small networks.

As we know the true ancestral states here, we can compare them to our
estimation. In this example, we see that the 95% prediction (ancestral state
reconstruction) intervals contain the true simulated value, at all ancestral nodes.

```@repl tree_trait
pred = predict(ancTrait1, interval=:prediction);
DataFrame(
  lower     = pred[1:7, :lower],                         # lower bound of 95% prediction interval
  trueValue = [3.312,4.438,3.922,3.342,2.564,1.315,2.0], # from sim1[:internalnodes] in next section
  upper     = pred[1:7, :upper]                          # upper bound
 )
```

### From estimated parameters

In real applications though, we do not have access to the true parameters of the
process that generated the data. We can estimate them using the previous function.
To fit a regular BM, we just need to do a regression of trait 1 against a simple
intercept:
```@example tree_trait
fitTrait1 = phylolm(@formula(trait1 ~ 1), dat, truenet)
nothing # hide
```

We can then apply the [`ancestralreconstruction`](@ref) function directly
to the fitted object:
```@example tree_trait
ancTrait1Approx = ancestralreconstruction(fitTrait1)
nothing # hide
```
The prediction intervals ignore the fact that we estimated the process
parameters, so they are less accurate and the function throws a warning.
The output is an object of the same [`ReconstructedStates`](@ref) type as earlier,
and the same extractors can be applied to it:
```@example tree_trait
R"svg"(name("ancestral1.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel = predict(ancTrait1Approx, text = true));
R"dev.off()" # hide
nothing # hide
```
![ancestral1](../assets/figures/ancestral1.svg)

For convenience, the two steps described above (fitting against the
intercept, and then do ancestral state reconstruction) can be done all at once
with a single call of the function [`ancestralreconstruction`](@ref) on a
DataFrame with the trait to reconstruct, and the tip labels:
```@example tree_trait
datTrait1 = DataFrame(trait1 = dat[:,:trait1], tipnames = dat[:,:tipnames])
ancTrait1Approx = ancestralreconstruction(datTrait1, truenet)
nothing # hide
```
```@example tree_trait
ancInt = predict(ancTrait1Approx, interval=:prediction, level=0.9, text=true, combine=true) 
R"svg"(name("ancestral2.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel = ancInt[!,[:nodenumber,:interval]]);
R"dev.off()" # hide
nothing # hide
```
![ancestral2](../assets/figures/ancestral2.svg)

This produces the exact same results. Here, we chose a `level` of 90% for the
plotted prediction intervals.

### Data imputation

Note that there is no theoretical difference between an internal node, for which
we could not measure the value of the trait, and a missing value at a tip of the
network. Consequently, the previous [`ancestralreconstruction`](@ref)
function can be used to do data imputation. To see this, let's add some missing
values in trait 1.
```@example tree_trait
allowmissing!(datTrait1, :trait1)
datTrait1[2, :trait1] = missing; # second row: for taxon C
ancTrait1Approx = ancestralreconstruction(datTrait1, truenet)
nothing # hide
```
```@example tree_trait
ancInt = predict(ancTrait1Approx, interval=:prediction, text=true, combine=true) 
R"svg"(name("ancestral3.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel = ancInt[!,[:nodenumber,:interval]]);
R"dev.off()" # hide
nothing # hide
```
![ancestral3](../assets/figures/ancestral3.svg)

A prediction interval is shown for the missing values.

### With known predictors

At this point, it might be tempting to apply this function to trait 3 as a
linear combination of trait 1 and a phylogenetic noise, to get a better
ancestral state reconstruction via using its correlation with trait 1.
However, this cannot be done directly:
```julia
ancTrait3 = ancestralreconstruction(fitTrait3) # throws an error
```
This is because the model to fit the trait (a regression with one
predictor and an intercept) used a predictor for which we don't know the
ancestral states. The regression model accounted for that.

The only option we have is to provide the function with the predictor's
ancestral states, if they are known. They are actually known in this
toy example because we generated the data ourselves (see next section),
so we can reconstruct our trait doing the following:
```@example tree_trait
ancTrait3 = ancestralreconstruction(fitTrait3, # model with estimated coefs etc.
  hcat(ones(7,1),
  [ 3.312, 4.438,  3.922,  3.342,  2.564,  1.315,  2.0], # from sim1[:internalnodes]
  [-3.62, -0.746, -2.217, -3.612, -2.052, -2.871, -2.0]) # from sim2[:internalnodes]
)
nothing # hide
```
```@example tree_trait
ancInt = predict(ancTrait3, interval=:prediction, text=true, combine=true) 
R"svg"(name("ancestral4.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, nodelabel = ancInt[!,[:nodenumber,:interval]]);
R"dev.off()" # hide
nothing # hide
```
![ancestral4](../assets/figures/ancestral4.svg)

where we provided the ancestral predictors as a matrix, containing the
intercept, and the known predictor at internal nodes. We must be very careful
with this function, as no check is done for the order of the predictors, or
the order of their values that must be the same as the internal nodes of the
phylogeny. As ancestral predictors are often unknown, the use of this
functionality is discouraged.

## Phylogenetic ANOVA

The [`phylolm`](@ref) function is based on the `lm` function
from [GLM](https://github.com/JuliaStats/GLM.jl). This means that it
inherits from most of its features, and in particular, it can handle formulas
with factors (discrete predictors) and interactions.
For example, in lizards, we might want to do a regression of toe length on
body length and the region where each species is found, where this region is coded
into 4 categories (say). We might also want to include an interaction effect
between body length and region.
(This model is just meant to show the possibilities of the function).

To illustrate the use of categorical predictors of particular interest in a
network with reticulations, let's assume that some transgressive evolution took
place after the hybridization event, so that species "A" and "B" have a larger
mean compared to the others
(see [^B18] for transgressive evolution after a reticulation event).
```@example tree_trait
δ = 5.0; # value of heterosis
underHyb = [n == "A" || n == "B" for n in dat[:,"tipnames"]] # tips under hybrid
underHyb
for i in 1:nrow(dat)
    underHyb[i] && (dat[i,:trait3] += δ) # add delta to tips A and B
end
nothing # hide
```
```@repl tree_trait
select(dat, [:trait3, :tipnames]) # trait3 changed: +5 added to A and B by previous loop
```
The categorical variable `underHyb` separates tips "A" and "B" from the others.
We need to consider it as a factor, not a numerical variable.
One way is to make it a vector of strings, as done below.
An alternative way would be to add and use the `CategoricalArrays` package,
then transform the column `underHyb` to be `categorical` (shown in commments).
```@example tree_trait
dat.underHyb = string.(underHyb); # adds a new column
# using CategoricalArrays
# transform!(dat, :underHyb => categorical, renamecols=false)
nothing # hide
```
```@repl tree_trait
dat
```
Now we can include this reticulation variable in the regression.
```@example tree_trait
fitTrait = phylolm(@formula(trait3 ~ trait1 + underHyb), dat, truenet)
```
In this case, the categorical variable indicating which tips are descendants
of the reticulation event is indeed relevant, and the transgressive evolution effect
is recovered.

This is a very simple example of how to include transgressive evolution,
but some general
functions to test for it, on networks with more than one hybrid, are also
available.

## Pagel's Lambda

One classical question about trait evolution is the amount of
"phylogenetic signal" in a trait or in the residuals of a linear relationship,
that is, the importance of the tree
structure to explain variation in the observed traits (or in the residuals).
One way of measuring that is to use
Pagel's lambda transformation of the branch lengths [^P99].
This model assumes a
BM on a tree where the internal branches are multiplied by a factor λ,
while the external branches are modified so that the total height of the tree is
constant. Hence, λ varies between 0 (the tree has no influence on
the data) and 1 (the tree is unchanged).
Using the same branch length transformations, this model can
be straightforwardly extended to phylogenetic networks.

This transformation assumes a time-consistent and ultrametric phylogeny,
in which all paths from the root to any tip has the same length: the "height"
of the phylogeny referred to above.

We can illustrate this with the predictor trait we used earlier. We use the
same function as before, only indicating the model we want to use:
```@example tree_trait
fitPagel = phylolm(@formula(trait1 ~ 1), dat, truenet, model="lambda")
```
As this trait 1 was indeed generated according to a plain BM on the phylogeny
(see next section),
the estimated λ should be close to 1. It can be extracted with function
`lambda_estim`:
```@repl tree_trait
lambda_estim(fitPagel)
```

For models in which the covariance is estimated, like Pagel's lambda,
model comparisons should use a likelihood ratio test with the function `lrtest`,
because the f-test (see below) is not applicable.

If the models being compared have different predictors, then models
should be fit with maximum likelihood instead of the default REML criterion
in order to do a likelihood ratio test: use option `reml=false` for this.

## Test of transgressive evolution

In the ANOVA section above, we showed how to include transgressive evolution
in a simple case, and we did so manually.
In general, transgressive evolution can be seen as a particular example
of a *shifted BM* on the phylogenetic network, in which the trait evolves
as a BM on the network, but *shifts* at a reticulation. The value of the shift
may be the same across all reticulations, or may differ between reticulations.

For identifiability reasons, each transgressive shift is applied to the
edge below a reticulation. In our network above, there is a single reticulation
and the edge below it is edge 6:

```@example tree_trait
R"svg"(name("truenet_with_numbers.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(truenet, useedgelength=true, showedgenumber=true);
R"dev.off()" # hide
nothing # hide
```
![truenet_with_numbers](../assets/figures/truenet_with_numbers.svg)

Let's assume we measured a trait that we hypothesized underwent a shift at
some or all ancestral reticulations. To test this hypothesis, we can use the 
custom columns of the [`PhyloNetworks.descendencematrix`](@extref), that can be
directly defined thanks to function [`descendencedataframe`](@ref).
```@repl tree_trait
df_shift = descendencedataframe(truenet) # regressors matching Hybrid Shifts
```
This creates a dataframe, with one column for each hybrid node
in the network, named according to the number of the edge after the
hybrid. In column `shift_6`, a row has a 0 if the corresponding species
is *not* a descendant of the reticulation, otherwise has the proportion of
its genome that was inherited from this reticulation. Here, A and B's ancestry
if fully inherited from edge 6, below the one reticulation in the network.

We can use the columns in this dataframe as regressors (predictors) in the
`phylolm` function. Their coefficients will measure the shift after each
reticulation.
In the example below, the species names are listed in a different order than in `df_shift`, and contained in a column called "species", to show how this is
handled to merge and then fit the data.

```@repl tree_trait
dat = DataFrame(  # trait data
  trait = [3.510, 2.195, 1.869, 4.839, 5.027, -0.679],
  species = ["O", "D", "C", "A", "B", "E"]);
dat = innerjoin(dat, df_shift, on = [:species => :tipnames]) # trait + shift predictors
fit_sh = phylolm(@formula(trait ~ shift_6), dat, truenet, tipnames=:species) # fit
```
Here, because there is only one hybrid in the network, we can directly
see whether the ancestral transgressive evolution is significant or not thanks to the
Student T-test on the coefficient associated with `shift_6`. In more
complex cases, it is possible to do a Fisher F-test, thanks to the `GLM`
function `ftest`.
```@example tree_trait
fit_null = phylolm(@formula(trait ~ 1), dat, truenet, tipnames=:species) # null (no shift)
ftest(fit_null, fit_sh)  # nested models
```
Here, this test is equivalent to the Fisher F test, and gives the same p-value.

!!! note "Warnings from GLM"
    A warning may appear, saying
    "Starting from GLM.jl 1.8, null model is defined as having no predictor at all when a model without an intercept is passed."
    - Why? `ftest` is inherited from the GLM package, which does not know that
      the intercept term is not a column of ones after transformation to remove
      the phylogenetic correlation. This is why `ftest` sends a warning for
      each model, when multiple models are compared.
    - These specific warnings can be ignored:
      * F values and p-values are correct
      * R² values are also correct: they are obtained with the
        `r2` function for phylogenetic linear models.
    A future version of the package will attempt to remove these warnings
    specifically.

Note that models need to be ordered by complexity, when given to `ftest`:
either from most complex to most simple, or from most simple to most complex.
In the output table, models are listed in the order in which they were given.
If the most complex model is given first, as done above, the table
lists the most complex H₁ (with shifts) first, and the null model H₀
is listed as the second model.

---

### References

[^B18]: Bastide, Solís-Lemus, Kriebel, Sparks, Ané (2018):
    Phylogenetic Comparative Methods for Phylogenetic Networks with Reticulations.
    Systematic Biology 67(5):800–820. doi:10.1093/sysbio/syy033

[^P99]: Pagel M (1999). Inferring the historical patterns of biological
    evolution. Nature. 401: 877–884. doi:10.1038/44766
