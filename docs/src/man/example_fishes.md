# Empirical example: Xiphophorus fishes

We reproduce here the analyses on *Xiphophorus* fishes from
[Bastide et al. (2018)](https://doi.org/10.1093/sysbio/syy033),
available on [dryad](https://doi.org/10.5061/dryad.nt2g6).

## Data loading

The data can be downloaded here:
- the time calibrated [network](https://datadryad.org/api/v2/files/343930/download);
- the [trait](https://datadryad.org/api/v2/files/343929/download)
  data on sword index and female preference for a sword,
  originally from [Cui et al. (2013)](https://doi.org/10.1111/evo.12099).

The files are also in the `examples` folder of the `PhyloTraits` package as
 `xiphophorus_networks_calibrated.tre` and `xiphophorus_morphology_Cui_etal_2013.csv`.

```@setup fish
using Logging
nowarninglogger = ConsoleLogger(stderr, Logging.Error)
mkpath("../assets/figures")
```

If not done already, load the packages needed for this analysis:
```@example fish
using PhyloNetworks, PhyloTraits, PhyloPlots
using RCall, CSV, DataFrames
using StatsAPI, StatsBase, StatsModels
name(x) = joinpath("..", "assets", "figures", x); # hide
nothing      # hide
```
then read in the networks data:
```@example fish
examples_path = joinpath(dirname(pathof(PhyloTraits)), "..", "examples");
topologies = readmultinewick(joinpath(examples_path, "xiphophorus_networks_calibrated.tre"));
net3 = topologies[3]; # we will use the network with 3 reticulations
R"svg"(name("net3.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net3; useedgelength=true); # topology + branch lengths
R"dev.off()" # hide
nothing      # hide
```

![net3](../assets/figures/net3.svg)

The file `examples/xiphophorus_networks_calibrated.tre` contains three
calibrated networks, with 0, 1 or 3 reticulations
(see [Bastide et al. (2018)](https://doi.org/10.1093/sysbio/syy033) for details).
In this tutorial, we use the network with 3 reticulations.

We can [`rotate!`](@extref PhyloNetworks.rotate!) some of the nodes
to avoid crossing edges and produce a better figure:

```@example fish
rotate!(net3, -4)
rotate!(net3, -5)
rotate!(net3, -6)
rotate!(net3, -14)
rotate!(net3, -16)
rotate!(net3, -17)
R"svg"(name("net3_rot.svg"), width=8, height=4) # hide
R"par"(mar=[0,0,0,0]) # hide
plot(net3; useedgelength=true); # topology + branch lengths
R"dev.off()" # hide
nothing      # hide
```

![net3](../assets/figures/net3_rot.svg)

We can then read in the trait data:
```@example fish
csvfile = joinpath(examples_path, "xiphophorus_morphology_Cui_etal_2013.csv")
dat = CSV.read(csvfile, DataFrame);
```

Sometimes the trait data and phylogeny have non-overlapping taxa sets.
`PhyloTraits` requires the data and phylogeny to have information on the same set of taxa.
We will delete the rows in the trait data for the taxa that are missing in the network.

```@repl fish
taxa_net = tiplabels(net3); # extract list of taxa
missing_rows = Integer[];
for i in reverse(1:nrow(dat))
    j = findfirst(isequal(dat.tipnames[i]), taxa_net)
    if isnothing(j) # taxon not found in network
        println("taxon ", dat.tipnames[i], " (row $i) not found in network")
        push!(missing_rows,i)
    end
end
dat = dat[setdiff(1:nrow(dat),missing_rows),:];
```
The snippet above should work fairly generically, assuming that the column with
the taxa is labelled `tipnames`.
Here, it tells us that
`taxon Xnezahualcoyotl (row 19) not found in network` and thus our code
(last line) removed it from our data frame `dat`.

## Ancestral state prediction

As [done in the tutorial](#From-estimated-parameters),
after fitting a model of trait evolution, we may wish to estimate the character
at various internal nodes to gain insight on the ancestral states.

For the ancestral states of the sword index:
```@example fish
dat_si = dat[:, [:tipnames,:sword_index]];
AS_si = ancestralreconstruction(dat_si, net3); # the ancestral state estimates
R"svg"(name("anc_fish_si_est.svg"), width=8, height=4) # hide
R"par(mar=c(.5,.5,.5,.5))"; # hide
plot(net3, nodelabel = predict(AS_si, text=true), xlim=[0,20]);
R"dev.off()" # hide
nothing      # hide
```
![anc_fish_si_est](../assets/figures/anc_fish_si_est.svg)

The prediction intervals ignore the fact that we estimated the process
parameters, so they are less accurate and the function throws a warning.
This warning will only show once in a given `julia` session.

To get the 90% credible interval in a text format to use as labels for plotting,
we can use option `text=true` in the `predict` function.
```@example fish
AS_si_int = predict(AS_si, interval=:prediction, level=0.90, text=true);
R"svg"(name("anc_fish_si_ci.svg"), width=8, height=4) # hide
R"par(mar=c(.5,.5,.5,.5))"; # hide
plot(net3, nodelabel=AS_si_int[!,[:nodenumber,:interval]], useedgelength=true);
R"dev.off()" # hide
nothing      # hide
```

![anc_fish_si_ci](../assets/figures/anc_fish_si_ci.svg)

Note that here, the first plot ignores branch lengths,
while the second one uses them (`useedgelength=true`).

As a consistency check, we can observe that the ancestral state prediction
interval falls within the extrema of the trait itself:
```@repl fish
AS_si_int[findfirst(AS_si_int.nodenumber .== -2),:interval]
extrema(dat[:,:sword_index])
```

We can do the same for the female preference for a sword:
```@example fish
dat_fp = dat[:, [:preference, :tipnames]];
with_logger(nowarninglogger) do # hide
global AS_fp=1 # hide
AS_fp = ancestralreconstruction(dat_fp, net3);
end # hide
R"svg"(name("anc_fish_fp_est.svg"), width=8, height=4) # hide
R"par(mar=c(.5,.5,.5,.5))"; # hide
plot(net3, nodelabel=predict(AS_fp,text=true), xlim=[0,20]);
R"dev.off()" # hide
nothing      # hide
```
![anc_fish_fp_est](../assets/figures/anc_fish_fp_est.svg)

```@example fish
AS_fp_int = predict(AS_fp,interval=:prediction,level=0.90,text=true);
R"svg"(name("anc_fish_fp_ci.svg"), width=8, height=4) # hide
R"par(mar=c(.5,.5,.5,.5))"; # hide
plot(net3, nodelabel=AS_fp_int[!,[:nodenumber,:interval]],
     useedgelength=true, xlim=[-3,26]);
R"dev.off()" # hide
nothing      # hide
```
![anc_fish_fp_ci](../assets/figures/anc_fish_fp_ci.svg)

```@repl fish
AS_fp_int[findfirst(AS_fp_int.nodenumber .== -2),:interval]
extrema(skipmissing(dat[:,:preference]))
```

## Phylogenetic signal: Pagel's lambda

We can use [Pagel's lambda](#Pagel's-Lambda) transformation to asses the
phylogenetic sigal.

```@repl fish
lambda_si = phylolm(@formula(sword_index ~ 1), dat, net3, model="lambda")
lambda_fp = phylolm(@formula(preference ~ 1),  dat, net3, model="lambda";
                    suppresswarnings=true)
```
On both traits we observe λ>1.0 when fitting the Pagel's lambda model;
consequently, we would interpret the patterns in trait data to have
high phylogenetic signal.
Note that here, we observe λ values greater than 1.0: the observed value is the
maximum value so that the transformation does not produce negative branch lengths.

!!! note "Warnings from PhyloTraits"
    When fitting a Pagel's lambda model on a network that is not time consistent,
    a warning will appear and notify the user that the node heights from the
    major tree within the network will be used for analysis.
    Pagel's lambda models assume time consistency, making it important to
    calibrate a network prior to model fitting.
    In our case, if we run `getnodeheights(net3)` we will get an error showing
    that the different path lengths that lead to the time inconsistency is
    relatively small. The small difference between paths likely results in a
    minor difference when using the major tree node heights instead of
    well-calibrated heights.
    ```@repl fish
    getnodeheights(net3)
    ```


## Phylogenetic regression of sword index versus preference

Phylogenetic regression can help us anwser the question:
does preference influence sword index?

```@example fish
R"svg"(name("sword_vs_preference.svg"), width=5, height=5) # hide
R"par"(mar=[3,2.5,.2,.2], las=1, mgp=[1.5,0.8,0]); # hide
R"plot"(dat.preference, dat.sword_index,
        xlab="female preference", ylab="sword index");
R"dev.off()" # hide
nothing      # hide
```
![sword_vs_preference](../assets/figures/sword_vs_preference.svg)

```@repl fish
fit_BM = phylolm(@formula(sword_index ~ preference), dat, net3)
fit_λ  = phylolm(@formula(sword_index ~ preference), dat, net3, model="lambda"; suppresswarnings=true)
lrtest(fit_BM,fit_λ)
```
On both Brownian motion and Pagel's lambda models,
we find a positive but statisitcally insignificant, relationship between mate preference and sword index.

Further, when comparing AIC values between the models, we may conclude that
including the extra λ parameter in the Pagel's lambda does not drastically
improve the model's ability to explain the patterns in the data.
Since the Brownian Motion model is nested within the Pagel's lambda model
(BM assumes λ=1.0), we can use a likelihood ratio test to more formally conclude
that the Pagel's lambda model does not significantly fit the data better than
the Brownain Motion alone.


## Transgressive evolution

To evaluate whether there was transgressive evolution that caused trait shifts
at reticulation events, we can fit and compare different models of
continuous trait evolution.

Here we compare three different models:
- `fit0`: no effect of reticulation events on the trait expectation;
- `fit1`: all three reticulation events lead to the same shift value on the trait expectation;
- `fit2`: each reticulation events has its own shift value on the trait expectation.

For the sword index, we get:
```@repl fish
df_shift = descendencedataframe(net3); # regressors matching Hybrid Shifts
dat3 = leftjoin(dat, df_shift, on = :tipnames); # add regressors to data
fit0 = phylolm(@formula(sword_index ~ 1),   dat3, net3) # no shift
fit1 = phylolm(@formula(sword_index ~ sum), dat3, net3) # same shift at hybrids
fit2 = phylolm(@formula(sword_index ~ shift_24 + shift_37 + shift_45),
               dat3, net3) # different shifts at hybrid nodes
with_logger(nowarninglogger) do # hide
ftest(fit0, fit1, fit2)
end # hide
```

From the p-values we might conclude that neither transgressive model can fit
the data significantly better than the model without trait shifts.

For the female preference, we get:
```@repl fish
fit0 = phylolm(@formula(preference ~ 1),   dat3, net3) # no shift
fit1 = phylolm(@formula(preference ~ sum), dat3, net3) # same shift at hybrids
fit2 = phylolm(@formula(preference ~ shift_24 + shift_37 + shift_45), dat3, net3) # different shifts at hybrid nodes
with_logger(nowarninglogger) do # hide
ftest(fit0, fit1, fit2)
ftest(fit0, fit2)
end # hide
```

The heterogeneous model gets some support, with effects in opposite directions
(some positive, some negative).
