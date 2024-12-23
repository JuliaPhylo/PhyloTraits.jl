###############################################################################

#= ------ roadmap of phylolm methods --------------

with or without within-species variation:
- phylolm(formula, dataframe, net; model="BM",...,withinspecies_var=false,...)
- phylolm(X,Y,net, model::ContinuousTraitEM; kwargs...)
  calls a function with or without within-species variation.

1. no measurement error in species means:
   - phylolm(model, X,Y,net, reml; kwargs...) dispatches based on model type
   - phylolm_lambda(X,Y,V,reml, gammas,times; ...)
   - phylolm_scalinghybrid(X,Y,net,reml, gammas; ...)

   helpers:
   - pgls(X,Y,V; ...) for vanilla BM, but called by others with fixed V_theta
   - logLik_lam(lambda, X,Y,V,gammas,times; ...)
   - logLik_lam_hyb(lambda, X,Y,net,gammas; ...)

2. with measurement error (within-species variation):
   - phylolm_wsp(model, X,Y,net, reml; kwargs...) dispatch based on model
     implemented for model <: BM only

   - phylolm_wsp(X,Y,V,reml, nonmissing,ind, counts,ySD, model_within)
   - phylolm_wsp(Xsp,Ysp,Vsp,reml, d_inv,RSS, n,p,a, model_within)
=#
"""
    phylolm(X::Matrix, Y::Vector, net::HybridNetwork, model::ContinuousTraitEM=BM(); kwargs...)

Return a [`PhyloNetworkLinearModel`](@ref) object.
This method is called by `phylolm(formula, data, network; kwargs...)`.
"""
function phylolm(
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    model::ContinuousTraitEM = BM();
    reml::Bool=true,
    nonmissing::BitArray{1}=trues(length(Y)),
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing,
    withinspecies_var::Bool=false,
    counts::Union{Nothing, Vector}=nothing,
    ySD::Union{Nothing, Vector}=nothing
)
    if withinspecies_var
        phylolm_wsp(model, X,Y,net, reml; nonmissing=nonmissing, ind=ind,
            ftolRel=ftolRel, xtolRel=xtolRel, ftolAbs=ftolAbs, xtolAbs=xtolAbs,
            counts=counts, ySD=ySD)
    else
        phylolm(model, X,Y,net, reml; nonmissing=nonmissing, ind=ind,
            ftolRel=ftolRel, xtolRel=xtolRel, ftolAbs=ftolAbs, xtolAbs=xtolAbs,
            startingValue=startingValue, fixedValue=fixedValue)
    end
end

function phylolm(
    ::BM,
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool;
    nonmissing::BitArray{1}=trues(length(Y)),
    ind::Vector{Int}=[0],
    kwargs...
)
    # BM variance covariance:
    # V_ij = expected shared time for independent genes in i & j
    V = sharedpathmatrix(net)
    linmod, Vy, RL, logdetVy = pgls(X,Y,V; nonmissing=nonmissing, ind=ind)
    return PhyloNetworkLinearModel(linmod, V, Vy, RL, Y, X,
                logdetVy, reml, ind, nonmissing, BM())
end

function phylolm(
    ::PagelLambda,
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool;
    nonmissing::BitArray{1}=trues(length(Y)),
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing
)
    # BM variance covariance
    V = sharedpathmatrix(net) # runs preorder!(net) by default
    gammas = getGammas(net)
    if istimeconsistent(net, false) # false: no need to preorder again
        times = getnodeheights(net, false)
    else
        @warn """The network is not time consistent (node heights are not well-defined).
        The network should be calibrated for this analysis, as the theory for Pagel's model
        assumes a time-consistent network.
        The analysis will use node heights based on the major tree, in the meantime.
        """
        times = getnodeheights_majortree(net, false; warn=false)
    end
    phylolm_lambda(X,Y,V,reml, gammas, times;
            nonmissing=nonmissing, ind=ind,
            ftolRel=ftolRel, xtolRel=xtolRel, ftolAbs=ftolAbs, xtolAbs=xtolAbs,
            startingValue=startingValue, fixedValue=fixedValue)
end


"""
    setGammas!(net, γ vector)

Set inheritance γ's of hybrid edges, using input vector for *major* edges.
Assume pre-order calculated already, with up-to-date field `vec_node`.
See [`getGammas`](@ref).

**Warning**: very different from [`PhyloNetworks.setgamma!`](@extref),
which focuses on a single hybrid event,
updates the field `ismajor` according to the new γ, and is not used here.

**Assumption**: each hybrid node has only 2 parents, a major and a minor parent
(according to the edges' field `ismajor`).
"""
function setGammas!(net::HybridNetwork, gammas::Vector)
    for (i,nod) in enumerate(net.vec_node)
        nod.hybrid || continue # skip tree nodes: nothing to do
        majorhyb = getparentedge(nod) # major
        minorhyb = getparentedgeminor(nod) # error if doesn't exit
        majorhyb.gamma = gammas[i]
        minorhyb.gamma = 1 - gammas[i]
    end
    return nothing
end

"""
    getGammas(net)

Vector of inheritance γ's of all major edges (tree edges and major hybrid edges),
ordered according to the pre-order index of their child node,
assuming this pre-order is already calculated
(with up-to-date field `vec_node`).
Here, a "major" edge is an edge with field `ismajor` set to true,
regardless of its actual γ (below, at or above 0.5).

See [`setGammas!`](@ref)
"""
function getGammas(net::HybridNetwork)
    gammas = ones(length(net.vec_node))
    for (i,node) in enumerate(net.vec_node)
        node.hybrid || continue # skip tree nodes: their gamma is already set to 1
        majorhybedge = getparentedge(node) # major
        gammas[i] = majorhybedge.gamma
    end
    return gammas
end

#= ScalingHybrid = BM but with optimized weights of hybrid edges:
minor edges have their original γ's changed to λγ. Same λ at all hybrids.
see Bastide (2017) dissertation, section 4.3.2 p.175, at
https://tel.archives-ouvertes.fr/tel-01629648
=#
function phylolm(
    ::ScalingHybrid,
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool;
    nonmissing::BitArray{1}=trues(length(Y)),
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing
)
    preorder!(net)
    gammas = getGammas(net)
    phylolm_scalinghybrid(X, Y, net, reml, gammas;
            nonmissing=nonmissing, ind=ind,
            ftolRel=ftolRel, xtolRel=xtolRel, ftolAbs=ftolAbs, xtolAbs=xtolAbs,
            startingValue=startingValue, fixedValue=fixedValue)
end

###############################################################################
## Fit BM

# Vanilla BM using covariance V. used for other models: V calculated beforehand
function pgls(
    X::Matrix,
    Y::Vector,
    V::MatrixTopologicalOrder;
    nonmissing::BitArray{1}=trues(length(Y)), # which tips are not missing?
    ind::Vector{Int}=[0]
)
    Vy = V[:tips] # extract tips matrix
    if (ind != [0]) Vy = Vy[ind, ind] end # re-order if necessary
    Vy = Vy[nonmissing, nonmissing]
    R = cholesky(Vy)
    RL = R.L
    # Fit with GLM.lm, and return quantities needed downstream
    return lm(RL\X, RL\Y), Vy, RL, logdet(R)
end

###############################################################################
## helper functions for lambda models

function maxLambda(times::Vector, V::MatrixTopologicalOrder)
    maskTips = indexin(V.tipnumbers, V.nodenumbers_toporder)
    maskNodes = indexin(V.internalnodenumbers, V.nodenumbers_toporder)
    return minimum(times[maskTips]) / maximum(times[maskNodes])
    # res = minimum(times[maskTips]) / maximum(times[maskNodes])
    # res = res * (1 - 1/5/maximum(times[maskTips]))
end

function transform_matrix_lambda!(
    V::MatrixTopologicalOrder,
    lam::AbstractFloat,
    gammas::Vector,
    times::Vector
)
    V.V .*= lam
    maskTips = indexin(V.tipnumbers, V.nodenumbers_toporder)
    for i in maskTips
        V.V[i, i] += (1-lam) * (gammas[i]^2 + (1-gammas[i])^2) * times[i]
    end
    #   V_diag = Matrix(Diagonal(diag(V.V)))
    #   V.V = lam * V.V .+ (1 - lam) .* V_diag
end

function logLik_lam(
    lam::AbstractFloat,
    X::Matrix,
    Y::Vector,
    V::MatrixTopologicalOrder,
    reml::Bool,
    gammas::Vector,
    times::Vector;
    nonmissing::BitArray{1}=trues(length(Y)), # Which tips are not missing ?
    ind::Vector{Int}=[0]
)
    # Transform V according to lambda
    Vp = deepcopy(V)
    transform_matrix_lambda!(Vp, lam, gammas, times)
    # Fit and take likelihood
    linmod, Vy, RL, logdetVy = pgls(X,Y,Vp; nonmissing=nonmissing, ind=ind)
    n = (reml ? dof_residual(linmod) : nobs(linmod))
    res = n*log(deviance(linmod)) + logdetVy
    if reml res += logdet(linmod.pp.chol); end
    return res
end

function phylolm_lambda(
    X::Matrix,
    Y::Vector,
    V::MatrixTopologicalOrder,
    reml::Bool,
    gammas::Vector,
    times::Vector;
    nonmissing::BitArray{1}=trues(length(Y)), # Which tips are not missing ?
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing
)
    if ismissing(fixedValue)
        # Find Best lambda using optimize from package NLopt
        opt = NLopt.Opt(:LN_BOBYQA, 1)
        NLopt.ftol_rel!(opt, ftolRel) # relative criterion
        NLopt.ftol_abs!(opt, ftolAbs) # absolute critetion
        NLopt.xtol_rel!(opt, xtolRel) # criterion on parameter value changes
        NLopt.xtol_abs!(opt, xtolAbs) # criterion on parameter value changes
        NLopt.maxeval!(opt, 1000) # max number of iterations
        NLopt.lower_bounds!(opt, 1e-100) # Lower bound
        # Upper Bound
        up = maxLambda(times, V)
        up = up-up/1000
        NLopt.upper_bounds!(opt, up)
        @info "Maximum lambda value to maintain positive branch lengths: " * @sprintf("%.6g", up)
        count = 0
        function fun(x::Vector{Float64}, g::Vector{Float64})
            x = convert(AbstractFloat, x[1])
            res = logLik_lam(x, X,Y,V, reml, gammas, times; nonmissing=nonmissing, ind=ind)
            count =+ 1
            #println("f_$count: $(round(res, digits=5)), x: $(x)")
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
        # Best value dans result
        res_lam = xmin[1]
    else
        res_lam = fixedValue
    end
    transform_matrix_lambda!(V, res_lam, gammas, times)
    linmod, Vy, RL, logdetVy = pgls(X,Y,V; nonmissing=nonmissing, ind=ind)
    res = PhyloNetworkLinearModel(linmod, V, Vy, RL, Y, X,
                logdetVy, reml, ind, nonmissing, PagelLambda(res_lam))
    return res
end

###############################################################################
## Fit scaling hybrid

function matrix_scalinghybrid(net::HybridNetwork, lam::AbstractFloat,
                              gammas::Vector)
    setGammas!(net, 1.0 .- lam .* (1. .- gammas))
    V = sharedpathmatrix(net)
    setGammas!(net, gammas)
    return V
end

function logLik_lam_hyb(
    lam::AbstractFloat,
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool,
    gammas::Vector;
    nonmissing::BitArray{1}=trues(length(Y)), # Which tips are not missing ?
    ind::Vector{Int}=[0]
)
    # Transform V according to lambda
    V = matrix_scalinghybrid(net, lam, gammas)
    # Fit and take likelihood
    linmod, Vy, RL, logdetVy = pgls(X,Y,V; nonmissing=nonmissing, ind=ind)
    n = (reml ? dof_residual(linmod) : nobs(linmod))
    res = n*log(deviance(linmod)) + logdetVy
    if reml res += logdet(linmod.pp.chol); end
    return res
end

function phylolm_scalinghybrid(
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool,
    gammas::Vector;
    nonmissing::BitArray{1}=trues(length(Y)), # Which tips are not missing ?
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing
)
    if ismissing(fixedValue)
        # Find Best lambda using optimize from package NLopt
        opt = NLopt.Opt(:LN_BOBYQA, 1)
        NLopt.ftol_rel!(opt, ftolRel) # relative criterion
        NLopt.ftol_abs!(opt, ftolAbs) # absolute critetion
        NLopt.xtol_rel!(opt, xtolRel) # criterion on parameter value changes
        NLopt.xtol_abs!(opt, xtolAbs) # criterion on parameter value changes
        NLopt.maxeval!(opt, 1000) # max number of iterations
        #NLopt.lower_bounds!(opt, 1e-100) # Lower bound
        #NLopt.upper_bounds!(opt, 1.0)
        count = 0
        function fun(x::Vector{Float64}, g::Vector{Float64})
            x = convert(AbstractFloat, x[1])
            res = logLik_lam_hyb(x, X, Y, net, reml, gammas; nonmissing=nonmissing, ind=ind)
            #count =+ 1
            #println("f_$count: $(round(res, digits=5)), x: $(x)")
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
        # Best value dans result
        res_lam = xmin[1]
    else
        res_lam = fixedValue
    end
    V = matrix_scalinghybrid(net, res_lam, gammas)
    linmod, Vy, RL, logdetVy = pgls(X,Y,V; nonmissing=nonmissing, ind=ind)
    res = PhyloNetworkLinearModel(linmod, V, Vy, RL, Y, X,
                logdetVy, reml, ind, nonmissing, ScalingHybrid(res_lam))
    return res
end

"""
    phylolm(f::StatsModels.FormulaTerm, fr::AbstractDataFrame, net::HybridNetwork; kwargs...)

Fit a phylogenetic linear regression model to data.
Return an object of type [`PhyloNetworkLinearModel`](@ref).
It contains a linear model from the GLM package, in `object.lm`, of type
[GLM.LinearModel](https://juliastats.org/GLM.jl/stable/api/#GLM.LinearModel).

## Arguments

* `f`: formula to use for the regression, see [StatsModels](https://juliastats.org/StatsModels.jl/stable/)
* `fr`: DataFrame containing the response values, predictor values, species/tip labels for each observation/row.
* `net`: phylogenetic network to use. Should have labelled tips.

Keyword arguments

* `model="BM"`: model for trait evolution (as a string)
  "lambda" (Pagel's lambda), "scalinghybrid" are other possible values
  (see [`ContinuousTraitEM`](@ref))
* `tipnames=:tipnames`: column name for species/tip-labels, represented
  as a symbol. For example, if the column containing the species/tip labels in
  `fr` is named "Species", then do `tipnames=:Species`.
* `no_names=false`: If `true`, force the function to ignore the tips names.
  The data is then assumed to be in the same order as the tips of the network.
  Default is false, setting it to true is dangerous, and strongly discouraged.
* `reml=true`: if `true`, use REML criterion ("restricted maximum likelihood")
  for estimating variance components, else use ML criterion.

The following tolerance parameters control the optimization of lambda if
`model="lambda"` or `model="scalinghybrid"`, and control the optimization of the
variance components if `model="BM"` and `withinspecies_var=true`.
* `fTolRel=1e-10`: relative tolerance on the likelihood value
* `fTolAbs=1e-10`: absolute tolerance on the likelihood value
* `xTolRel=1e-10`: relative tolerance on the parameter value
* `xTolAbs=1e-10`: absolute tolerance on the parameter value

* `startingValue=0.5`: If `model`="lambda" or "scalinghybrid", this
  provides the starting value for the optimization in lambda.
* `fixedValue=missing`: If `model`="lambda" or "scalinghybrid", and
  `fixedValue` is a number, then lambda is set to this number and is not optimized.
* `withinspecies_var=false`: If `true`, fits a within-species variation model.
  Currently only implemented for `model`="BM".
* `y_mean_std::Bool=false`: If `true`, and `withinspecies_var=true`, then accounts for
  within-species variation, using species-level statistics provided in `fr`.

## Methods applied to fitted models

To access the response values, do `response(object)`.
To access the model matrix, do `modelmatrix(object)`.
To access the model formula, do `formula(object)`.

## Within-species variation

For a high-level description, see [`PhyloNetworkLinearModel`](@ref).
To fit a model with within-species variation in the response variable,
either of the following must be provided in the data frame `fr`:

(1) Individual-level data: There should be columns for response, predictors, and
species/tip-labels. Every row should correspond to an individual observation.
At least one species must be represented by two or more individuals.

(2) Species-level statistics: There should be columns for mean response, predictors,
species/tip-labels, species sample-sizes (number of individuals for each species),
and species standard deviations (standard deviations of the response values
by species). Every row should correspond to a species: each species should be
represented by a unique row. The column names for species sample-sizes and
species standard deviations are expected to be "[response column name]\\_n"
and "[response column name]\\_sd". For example, if the response column name is "y",
then the column names should be "y\\_n" and "y\\_sd" for the sample-sizes and
standard deviations.

Regardless of whether the data provided follows (1) or (2),
`withinspecies_var` should be set to true.
If the data provided follows (2), then `y_mean_std` should be set to false.

## Within-species variation in predictors

The model assumes *no* within-species variation in predictors, because it aims to
capture the evolutionary (historical, phylogenetic) relationship between the
predictors and the response, not the within-species (present-day, or phenotypic)
relationship.

If a within-species variation model is fitted on individual-level data, and
if there are individuals within the same species with different values for
the same predictor, these values are all replaced by the mean predictor value
for all the individuals in that species.
For example, suppose there are 3 individuals in a given species, and that their
predictor values are (x₁=3, x₂=6), (x₁=4, x₂=8) and (x₁=2, x₂=1). Then the predictor
values for these 3 individuals are each replaced by (x₁=(3+4+2)/3, x₂=(6+8+1)/3)
before model fitting. If a fourth individual had data (x₁=10, x₂=missing), then that
individual would be ignored for any model using x₂, and would not contribute any
information to its species data for these models.

## Missing data

Rows with missing data for either the response or the predictors are omitted from
the model-fitting. There should minimally be columns for response, predictors,
species/tip-labels. As detailed above, additional columns may be required for fitting
within-species variation. Missing data in the columns for species names,
species standard deviation / sample sizes (if used) will throw an error.

## See also

[`rand`](@ref), [`ancestralreconstruction`](@ref), [`PhyloNetworks.vcv`](@extref)

## Examples: Without within-species variation

We first load data from the package and fit the default BM model.

```jldoctest phylolmdoc
julia> phy = readnewick(joinpath(dirname(pathof(PhyloTraits)), "..", "examples", "caudata_tree.txt"));

julia> using DataFrames, CSV # to read data file, next

julia> dat = CSV.read(joinpath(dirname(pathof(PhyloTraits)), "..", "examples", "caudata_trait.txt"), DataFrame);

julia> using StatsModels # for stat model formulas

julia> fitBM = phylolm(@formula(trait ~ 1), dat, phy; reml=false);

julia> fitBM # Shows a summary
PhyloNetworkLinearModel

Formula: trait ~ 1

Model: Brownian motion

Parameter Estimates, using ML:
phylogenetic variance rate: 0.00294521

Coefficients:
─────────────────────────────────────────────────────────────────────
             Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
─────────────────────────────────────────────────────────────────────
(Intercept)  4.679    0.330627  14.15    <1e-31    4.02696    5.33104
─────────────────────────────────────────────────────────────────────
Log Likelihood: -78.9611507833
AIC: 161.9223015666
```

We can extra parameters, likelihood, AIC etc.
```jldoctest phylolmdoc
julia> round(sigma2_phylo(fitBM), digits=6) # rounding for jldoctest convenience
0.002945

julia> round(mu_phylo(fitBM), digits=4)
4.679

julia> using StatsBase # for aic() stderror() loglikelihood() etc.

julia> round(loglikelihood(fitBM), digits=10)
-78.9611507833

julia> round(aic(fitBM), digits=10)
161.9223015666

julia> round(aicc(fitBM), digits=10)
161.9841572367

julia> round(bic(fitBM), digits=10)
168.4887090241

julia> round.(coef(fitBM), digits=4)
1-element Vector{Float64}:
 4.679

julia> confint(fitBM) # 95% (default) confidence interval for the coefficient(s)
1×2 Matrix{Float64}:
 4.02696  5.33104

julia> abs(round(r2(fitBM), digits=10)) # absolute value for jldoctest convenience
0.0

julia> abs(round(adjr2(fitBM), digits=10))
0.0

julia> round.(vcov(fitBM), digits=6) # variance-covariance of estimated parameters: squared standard error
1×1 Matrix{Float64}:
 0.109314
```

The residuals are the variance not explained by predictors.
The phylogenetic correlation modelled by the BM is about them.
The trait may have 2 sources of phylogenetic signal: from the predictor with
which it the response may be associated, and from the residuals.

```jldoctest phylolmdoc
julia> round.(residuals(fitBM), digits=6)
197-element Vector{Float64}:
 -0.237648
 -0.357937
 -0.159387
 -0.691868
 -0.323977
 -0.270452
 -0.673486
 -0.584654
 -0.279882
 -0.302175
  ⋮
 -0.777026
 -0.385121
 -0.443444
 -0.327303
 -0.525953
 -0.673486
 -0.603158
 -0.211712
 -0.439833

julia> round.(response(fitBM), digits=5)
197-element Vector{Float64}:
 4.44135
 4.32106
 4.51961
 3.98713
 4.35502
 4.40855
 4.00551
 4.09434
 4.39912
 4.37682
 ⋮
 3.90197
 4.29388
 4.23555
 4.3517
 4.15305
 4.00551
 4.07584
 4.46729
 4.23917

julia> round.(predict(fitBM), digits=5)
197-element Vector{Float64}:
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 ⋮
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
 4.679
```

## Examples: With within-species variation (two different input formats shown)

We use a smaller network here.
We can input data as 1 row per individual, multiple rows per species:

```jldoctest phylolmdoc
julia> net = readnewick("((((D:0.4,C:0.4):4.8,((A:0.8,B:0.8):2.2)#H1:2.2::0.7):4.0,(#H1:0::0.3,E:3.0):6.2):2.0,O:11.2);");

julia> df = DataFrame( # individual-level observations
           species = repeat(["D","C","A","B","E","O"],inner=3),
           trait1 = [4.08298,4.08298,4.08298,3.10782,3.10782,3.10782,2.17078,2.17078,2.17078,1.87333,1.87333,
              1.87333,2.8445,2.8445,2.8445,5.88204,5.88204,5.88204],
           trait2 = [-7.34186,-7.34186,-7.34186,-7.45085,-7.45085,-7.45085,-3.32538,-3.32538,-3.32538,-4.26472,
              -4.26472,-4.26472,-5.96857,-5.96857,-5.96857,-1.99388,-1.99388,-1.99388],
           trait3 = [18.8101,18.934,18.9438,17.0687,17.0639,17.0732,14.4818,14.1112,14.2817,13.0842,12.9562,
              12.9019,15.4373,15.4075,15.4317,24.2249,24.1449,24.1302]);

julia> m1 = phylolm(@formula(trait3 ~ trait1), df, net;
                    tipnames=:species, withinspecies_var=true)
PhyloNetworkLinearModel

Formula: trait3 ~ 1 + trait1

Model: Brownian motion

Parameter Estimates, using REML:
phylogenetic variance rate: 0.156188
within-species variance: 0.0086343

Coefficients:
──────────────────────────────────────────────────────────────────────
               Coef.  Std. Error     t  Pr(>|t|)  Lower 95%  Upper 95%
──────────────────────────────────────────────────────────────────────
(Intercept)  9.65347    1.3066    7.39    0.0018    6.02577   13.2812
trait1       2.30358    0.276163  8.34    0.0011    1.53683    3.07033
──────────────────────────────────────────────────────────────────────
Log Likelihood: 1.9446255188
AIC: 4.1107489623
```

Alternatively, we can input the data as 1 row per species and 2 extra columns:
standard deviation of the response trait among individuals in the same species,
and number of individuals per species for which we calculated this SD.
The result is the same.

```jldoctest phylolmdoc
julia> df_r = DataFrame( # species-level statistics (sample means, standard deviations)
           species = ["D","C","A","B","E","O"],
           trait1 = [4.08298,3.10782,2.17078,1.87333,2.8445,5.88204],
           trait2 = [-7.34186,-7.45085,-3.32538,-4.26472,-5.96857,-1.99388],
           trait3 = [18.896,17.0686,14.2916,12.9808,15.4255,24.1667],
           trait3_sd = [0.074524,0.00465081,0.185497,0.0936,0.0158379,0.0509643],
           trait3_n = [3, 3, 3, 3, 3, 3]);

julia> m2 = phylolm(@formula(trait3 ~ trait1), df_r, net;
                tipnames=:species, withinspecies_var=true, y_mean_std=true)
PhyloNetworkLinearModel

Formula: trait3 ~ 1 + trait1

Model: Brownian motion

Parameter Estimates, using REML:
phylogenetic variance rate: 0.15618
within-species variance: 0.0086343

Coefficients:
──────────────────────────────────────────────────────────────────────
               Coef.  Std. Error     t  Pr(>|t|)  Lower 95%  Upper 95%
──────────────────────────────────────────────────────────────────────
(Intercept)  9.65342    1.30657   7.39    0.0018    6.02582   13.281
trait1       2.30359    0.276156  8.34    0.0011    1.53686    3.07032
──────────────────────────────────────────────────────────────────────
Log Likelihood: 1.9447243714
AIC: 4.1105512573
```
"""
function phylolm(
    f::StatsModels.FormulaTerm,
    fr::AbstractDataFrame,
    net::HybridNetwork;
    model::AbstractString="BM",
    tipnames::Symbol=:tipnames,
    no_names::Bool=false,
    reml::Bool=true,
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    startingValue::Real=0.5,
    fixedValue::Union{Real,Missing}=missing,
    withinspecies_var::Bool=false,
    y_mean_std::Bool=false
)
    # Match the tips names: make sure that the data provided by the user will
    # be in the same order as the ordered tips in matrix V.
    preorder!(net)
    if no_names # The names should not be taken into account.
        ind = [0]
        @info """As requested (no_names=true), I am ignoring the tips names
             in the network and in the dataframe."""
    else
        nodatanames = !any(DataFrames.propertynames(fr) .== tipnames)
        nodatanames && any(tiplabels(net) == "") &&
            error("""The network provided has no tip names, and the input dataframe has no
                  column "$tipnames" for species names, so I can't match the data on the network
                  unambiguously. If you are sure that the tips of the network are in the
                  same order as the values of the dataframe provided, then please re-run
                  this function with argument no_name=true.""")
        any(tiplabels(net) == "") &&
            error("""The network provided has no tip names, so I can't match the data
                  on the network unambiguously. If you are sure that the tips of the
                  network are in the same order as the values of the dataframe provided,
                  then please re-run this function with argument no_name=true.""")
        nodatanames &&
            error("""The input dataframe has no column "$tipnames" for species names, so I can't
                  match the data on the network unambiguously. If you are sure that the
                  tips of the network are in the same order as the values of the dataframe
                  provided, then please re-run this function with argument no_name=true.""")
        ind = indexin(fr[!, tipnames], tiplabels(net))
        any(isnothing, ind) &&
            error("""Tips with data are not in the network: $(fr[isnothing.(ind), tipnames])
                  please provide a larger network including these tips.""")
        ind = convert(Vector{Int}, ind) # Int, not Union{Nothing, Int}
        if length(unique(ind)) == length(ind)
            withinspecies_var && !y_mean_std &&
            error("""for within-species variation, at least 1 species must have at least 2 individuals.
                  did you mean to use option "y_mean_std=true" perhaps?""")
        else
            (!withinspecies_var || y_mean_std) &&
            error("""Some tips have data on multiple rows.""")
        end
    end
    # Find the regression matrix and response vector
    data, nonmissing = StatsModels.missing_omit(StatsModels.columntable(fr), f)
    sch = StatsModels.schema(f, data)
    f = StatsModels.apply_schema(f, sch, PhyloNetworkLinearModel)
    mf = ModelFrame(f, sch, data, PhyloNetworkLinearModel)
    mm = StatsModels.ModelMatrix(mf)
    Y = StatsModels.response(mf)

    if withinspecies_var && y_mean_std
        # find columns in data frame for: # of individuals from each species
        colname = Symbol(String(mf.f.lhs.sym)*"_n")
        any(DataFrames.propertynames(fr) .== colname) ||
            error("expected # of individuals in column $colname, but no such column was found")
        counts  = fr[nonmissing,colname]
        all(!ismissing, counts) || error("some num_individuals values are missing, column $colname")
        all(x -> x>0, counts) || error("some species have 0 or <0 num_individuals, column $colname")
        all(isfinite.(counts))|| error("some species have infinite num_individuals, column $colname")
        # find sample SDs corresponding to the response mean in each species
        colname = Symbol(String(mf.f.lhs.sym)*"_sd")
        any(DataFrames.propertynames(fr) .== colname) ||
            error("expected the response's SD (per species) in column $colname, but no such column was found")
        ySD = fr[nonmissing,colname]
        all(!ismissing, ySD) || error("some SD values are missing, column $colname")
        all(x -> x≥0, ySD) || error("some SD values are negative, column $colname")
        all(isfinite.(ySD))|| error("some SD values are infinite, column $colname")
    else
        counts = nothing
        ySD = nothing
    end

    withinspecies_var && model != "BM" &&
        error("within-species variation is not implemented for non-BM models")
    modeldic = Dict("BM" => BM(),
                    "lambda" => PagelLambda(),
                    "scalinghybrid" => ScalingHybrid())
    haskey(modeldic, model) || error("phylolm is not defined for model $model.")
    modelobj = modeldic[model]

    res = phylolm(mm.m, Y, net, modelobj; reml=reml, nonmissing=nonmissing, ind=ind,
                  ftolRel=ftolRel, xtolRel=xtolRel, ftolAbs=ftolAbs, xtolAbs=xtolAbs,
                  startingValue=startingValue, fixedValue=fixedValue,
                  withinspecies_var=withinspecies_var, counts=counts, ySD=ySD)
    res.formula = f
    return res
end

###############################################################################
#  within-species variation (including measurement error)
###############################################################################

function phylolm_wsp(
    ::BM,
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool;
    nonmissing::BitArray{1}=trues(length(Y)), # which individuals have non-missing data?
    ind::Vector{Int}=[0],
    ftolRel::AbstractFloat=fRelTr,
    xtolRel::AbstractFloat=xRelTr,
    ftolAbs::AbstractFloat=fAbsTr,
    xtolAbs::AbstractFloat=xAbsTr,
    counts::Union{Nothing, Vector}=nothing,
    ySD::Union{Nothing, Vector}=nothing
)
    V = sharedpathmatrix(net)
    phylolm_wsp(X,Y,V, reml, nonmissing,ind,
        ftolRel, xtolRel, ftolAbs, xtolAbs,
        counts,ySD)
end

#= notes about missing data: after X and Y produced by stat formula:
- individuals with missing data (in response or any predictor)
  already removed from both X and Y
- V has all species: some not listed, some listed but without any data
- nonmissing and ind correspond to the original rows in the data frame,
  including those with some missing data, so:
  * nonmissing has length >= length of Y
  * sum(nonmissing) = length of Y
- V[:tips][ind,ind][nonmissing,nonmissing] correspond to the data rows

extra problems:
- a species may be listed 1+ times in ind, but not in ind[nonmissing]
- ind and nonmissing need to be converted to the species level, alongside Y
=#
function phylolm_wsp(
    X::Matrix,
    Y::Vector,
    V::MatrixTopologicalOrder,
    reml::Bool,
    nonmissing::BitArray{1},
    ind::Vector{Int},
    ftolRel::AbstractFloat,
    xtolRel::AbstractFloat,
    ftolAbs::AbstractFloat,
    xtolAbs::AbstractFloat,
    counts::Union{Nothing, Vector},
    ySD::Union{Nothing, Vector}
)
    n_coef = size(X, 2) # no. of predictors
    individualdata = isnothing(counts)
    xor(individualdata, isnothing(ySD)) &&
        error("counts and ySD must be both nothing, or both vectors")
    if individualdata
        # get species means for Y and X, the within-species residual ss
        ind_nm = ind[nonmissing] # same length as Y
        ind_sp = unique(ind_nm)
        n_sp = length(ind_sp) # number of species with data
        n_tot = length(Y)     # total number of individuals with data
        d_inv = zeros(n_sp)
        Ysp = Vector{Float64}(undef,n_sp) # species-level mean Y response
        Xsp = Matrix{Float64}(undef,n_sp,n_coef)
        RSS = 0.0  # residual sum-of-squares within-species
        for (i0,iV) in enumerate(ind_sp)
            iii = findall(isequal(iV), ind_nm)
            n_i = length(iii) # number of obs for species of index iV in V
            d_inv[i0] = 1/n_i
            Xsp[i0, :] = mean(X[iii, :], dims=1) # ideally, all have same Xs
            ymean = mean(Y[iii])
            Ysp[i0] = ymean
            RSS += sum((Y[iii] .- ymean).^2)
        end
        Vsp = V[:tips][ind_sp,ind_sp]
        # redefine "ind" and "nonmissing" at species level. ind = index of species
        # in tiplabels(net), in same order in which species come in means Ysp.
        # nonmissing: no need to list species with no data
        ind = ind_sp
        nonmissing = trues(n_sp)
    else # group means and sds for response variable were passed in
        n_sp = length(Y)
        n_tot = sum(counts)
        d_inv = 1.0 ./ counts
        Ysp = Y
        Xsp = X
        RSS = sum((ySD .^ 2) .* (counts .- 1.0))
        ind_nm = ind[nonmissing]
        Vsp = V[:tips][ind_nm,ind_nm]
    end

    model_within, RL = withinsp_varianceratio(Xsp,Ysp,Vsp, reml,
        ftolRel, xtolRel, ftolAbs, xtolAbs,
        d_inv, RSS, n_tot, n_coef, n_sp)
    η = model_within.optsum.final[1]
    Vm = Vsp + η * Diagonal(d_inv)
    m = PhyloNetworkLinearModel(lm(RL\Xsp, RL\Ysp), V, Vm, RL, Ysp, Xsp,
            2*logdet(RL), reml, ind, nonmissing, BM(), model_within)
    return m
end

# the method below takes in "clean" X,Y,V: species-level means, no missing data,
#     matching order of species in X,Y and V, no extra species in V.
# given V & η: analytical formula for σ² estimate
# numerical optimization of η = σ²within / σ²
function withinsp_varianceratio(
    X::Matrix,
    Y::Vector,
    V::Matrix,
    reml::Bool,
    ftolRel::AbstractFloat,
    xtolRel::AbstractFloat,
    ftolAbs::AbstractFloat,
    xtolAbs::AbstractFloat,
    d_inv::Vector,
    RSS::Float64,
    ntot::Real,
    ncoef::Int64,
    nsp::Int64,
    model_within::Union{Nothing, WithinSpeciesCTM}=nothing
)
    RL = cholesky(V).L
    lm_sp = lm(RL\X, RL\Y)
    if model_within === nothing
        # create model_within with good starting values
        s2start = GLM.dispersion(lm_sp, false) # sqr=false: deviance/dof_residual
        # this is the REML, not ML estimate, which would be deviance/nobs
        s2withinstart = RSS/(ntot-nsp)
        ηstart = s2withinstart / s2start
        optsum = OptSummary([ηstart], [1e-100], :LN_BOBYQA; initial_step=[0.01],
            ftol_rel=ftolRel, ftol_abs=ftolAbs, xtol_rel=xtolRel, xtol_abs=[xtolAbs])
        optsum.maxfeval = 1000
        model_within = WithinSpeciesCTM([s2withinstart], [s2start], d_inv, RSS, optsum)
    else
        optsum = model_within.optsum
        # fixit: I find this option dangerous (and not used). what if the
        # current optsum has 2 parameters instead of 1, or innapropriate bounds, etc.?
        # We could remove the option to provide a pre-built model_within
    end
    opt = Opt(optsum)
    Ndof = (reml ? ntot - ncoef : ntot )
    wdof = ntot - nsp
    Vm = similar(V) # scratch space for repeated usage
    function logliksigma(η) # returns: -2loglik, estimated sigma2, and more
        Vm .= V + η * Diagonal(d_inv)
        Vmchol = cholesky(Vm) # LL' = Vm
        RL = Vmchol.L
        lm_sp = lm(RL\X, RL\Y)
        σ² = (RSS/η + deviance(lm_sp))/Ndof
        # n2ll = -2 loglik except for Ndof*log(2pi) + sum log(di) + Ndof
        n2ll = Ndof * log(σ²) + wdof * log(η) + logdet(Vmchol)
        if reml
            n2ll += logdet(lm_sp.pp.chol) # log|X'Vm^{-1}X|
        end
        #= manual calculations without cholesky
        Q = X'*(Vm\X);  β = Q\(X'*(Vm\Ysp));  r = Y-X*β
        val =  Ndof*log(σ²) + ((RSS/η) + r'*(Vm\r))/σ² +
            (ntot-ncoef)*log(η) + logabsdet(Vm)[1] + logabsdet(Q)[1]
        =#
        return (n2ll, σ², Vmchol)
    end
    obj(x, g) = logliksigma(x[1])[1] # x = [η]
    NLopt.min_objective!(opt, obj)
    fmin, xmin, ret = NLopt.optimize(opt, optsum.initial)
    optsum.feval = opt.numevals
    optsum.final = xmin
    optsum.fmin = fmin
    optsum.returnvalue = ret
    # save the results
    η = xmin[1]
    (n2ll, σ², Vmchol) = logliksigma(η)
    model_within.wsp_var[1] = η*σ²
    model_within.bsp_var[1] = σ²
    return model_within, Vmchol.L
end

###############################################################################
#= Model comparisons

isnested: borrowed from GLM.issubmodel (as of commit 504e5186c87)
https://github.com/JuliaStats/GLM.jl/blob/master/src/ftest.jl#L11
To avoid comparing the coefnames and to be less restrictive, we compare the
design matrices. For example: Xsmall = [x1-x2 x1-x3] is nested in Xbig = [x1 x2 x3].
We check that there exists B such that Xsmall = Xbig * B, or rather, that
min_B norm(Xbig*B - Xsmall) ≈ 0 . For the math of this minimization problem,
see https://github.com/JuliaStats/GLM.jl/pull/197#issuecomment-331136617

When isnested() exists in GLM, check to see if we should improve further.
=#
"""
    isnested(m1::PhyloNetworkLinearModel, m2::PhyloNetworkLinearModel)
    isnested(m1::ContinuousTraitEM, m2::ContinuousTraitEM)

True if `m1` is nested in `m2`, false otherwise.
Models fitted with different criteria (ML and REML) are not nested.
Models with different predictors (fixed effects) must be fitted with ML to be
considered nested.
"""
function StatsModels.isnested(m1m::PhyloNetworkLinearModel, m2m::PhyloNetworkLinearModel; atol::Real=0.0)
    if !(nobs(m1m) ≈ nobs(m2m))
        @error "Models must have the same number of observations"
        return false
    end
    # exact same response? (before phylogenetic transformation)
    if m1m.Y != m2m.Y
        @error "Models must fit the same response"
        return false
    end
    # same criterion?
    if xor(m1m.reml, m2m.reml)
        @error "Models must be fitted with same criterion (both ML or both REML)"
        return false
    end
    # same within-species variation? e.g. same assumption of species means
    # this check should be useless, because same Y so same # species, and same
    # nobs so same # individuals. But 1 ind/species w/o within-species variation,
    # and at least 1 species with 2+ inds/species w/ within-species variation.
    xor(isnothing(m1m.model_within), isnothing(m2m.model_within)) && return false
    # nesting of fixed effects: is X1 = X2*B for some B?
    X1 = m1m.X
    np1 = size(X1, 2)
    X2 = m2m.X
    np2 = size(X2, 2)
    np1 > np2 && return false # if X1 has more predictors, can't be nested in X2
    # is mininum_B norm X2*B - X1 ≈ 0 ?
    rtol = Base.rtoldefault(eltype(X1))
    norm(view(qr(X2).Q' * X1, np2 + 1:size(X2,1), :)) < max(atol, rtol*norm(X1)) ||
        return false
    # ML (not REML) if different fixed effects
    sameFE = (X1 == X2) # exact equality okay here
    if !sameFE && (m1m.reml || m2m.reml)
        @error "Models should be fitted with ML to do a likelihood ratio test with different predictors"
        return false
    end
    # nesting of phylogenetic variance models
    return isnested(m1m.evomodel, m2m.evomodel)
end

isnested(::T,::T) where T <: ContinuousTraitEM = true
isnested(::BM,::Union{PagelLambda,ScalingHybrid}) = true
isnested(::Union{PagelLambda,ScalingHybrid}, ::BM) = false
isnested(::ScalingHybrid,::PagelLambda) = false
isnested(::PagelLambda,::ScalingHybrid) = false

#= ANOVA using ftest from GLM - need version 0.8.1
 As of GLM v1.8, ftest throws a warning on typical BM models, one per model:
 "Starting from GLM.jl 1.8, null model is defined as having no predictor at all when a model without an intercept is passed."
 This is because after transforming the data to de-correlate the residuals,
 the transformed intercept vector is not proportional to the constant vector 1.
 The warning is from: ftest → r2(phylomodel.lm) → nulldeviance(phylomodel.lm) → warning.
 R² by GLM is wrong: assume *no* intercept, and are based on the transformed data.
 R² corrected them below: r2(phylomodel) reimplemented here.
 But nulldeviance(phylomodel) does *not* call nulldeviance(phylomodel.lm),
 instead re-implemented here to use the intercept properly.
 Keep the warnings: unless they can be suppressed with specificity
 Ideally: modify `ftest` here or in GLM.
=#
function GLM.ftest(objs::PhyloNetworkLinearModel...; kwargs...)
    if !all( isa(o.evomodel,BM) && isnothing(o.model_within) for o in objs)
        throw(ArgumentError("""F test is only valid for the vanilla BM model.
        Use a likelihood ratio test instead with function `lrtest`."""))
    end
    objslm = [obj.lm for obj in objs]
    resGLM = ftest(objslm...; kwargs...)
    resGLM.r2 = r2.(objs)
    return resGLM
end
## ANOVA: old version - kept for tests purposes - do not export
"""
    anova(objs::PhyloNetworkLinearModel...)

Takes several nested fits of the same data, and computes the F statistic for each
pair of models.

The fits must be results of function [`phylolm`](@ref) called on the same
data, for models that have more and more effects.

Returns a DataFrame object with the anova table.
"""
function anova(objs::PhyloNetworkLinearModel...)
    anovaTable = Array{Any}(undef, length(objs)-1, 6)
    ## Compute binary statistics
    for i in 1:(length(objs) - 1)
      anovaTable[i, :] = anovaBin(objs[i], objs[i+1])
    end
    ## Transform into a DataFrame
    anovaTable = DataFrame(anovaTable,
        [:dof_res, :RSS, :dof, :SS, :F, Symbol("Pr(>F)")])
    return(anovaTable)
end

function anovaBin(obj1::PhyloNetworkLinearModel, obj2::PhyloNetworkLinearModel)
    length(coef(obj1)) < length(coef(obj2)) || error("Models must be nested, from the smallest to the largest.")
    ## residuals
    dof2 = dof_residual(obj2)
    dev2 = deviance(obj2, Val(true))
    ## reducted residuals
    dof1 = dof_residual(obj1) - dof2
    dev1 = deviance(obj1, Val(true)) - dev2
    ## Compute statistic
    F = (dev1 / dof1) / (dev2 / dof2)
    pval = GLM.ccdf.(GLM.FDist(dof1, dof2), F) # ccdf and FDist from Distributions, used by GLM
    return([dof2, dev2, dof1, dev1, F, pval])
end

###############################################################################
## Regression matrices for ANOVA / shifts

"""
    descendencedataframe(node::Vector{Node}, net::HybridNetwork; checkpreorder=true)
    descendencedataframe(edge::Vector{Edge}, net::HybridNetwork; checkpreorder=true)
    descendencedataframe(net::HybridNetwork; checkpreorder=true, which=:allhybrids)

Data frame containing the genomic proportion inherited by each taxon in `net`,
from each node or edge in the input vector, or from each hybrid node by default
in the third method. The data frame has 1 row per tip (taxon) in the network
and the following columns:
- 1 column per edge or node, with columns named according to the pattern
  shift_{edge_number}" where `edge_number` is the number of the edge associated
  with an input edge or node (in which case the parent edge is used)
- 1 additional column labeled `tipnames`, containing the tip labels.

The `shift_*` columns in this data frame can be used as regressor vectors
associated with shifts on input `edge`s or on edges that are above input `node`s.
With option `which=:allhybrids` in last method, each shift column is associated
with a hybrid node in `net`, to model a shift on the edge that is immediately
below the hybrid node. This can be used to test for transgressive evolution:
see examples below.

These methods use [`PhyloNetworks.descendencematrix`](@extref), so
`net` might be modified to store a vector of its nodes sorted in a pre-order.

# Examples

```jldoctest descendence
julia> net = readnewick("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);");

julia> preorder!(net)

julia> using PhyloPlots

julia> plot(net, shownodenumber=true); # to locate nodes

julia> nodes_shifts = indexin([1,-5], [n.number for n in net.node]) # Put a shift on edges ending at nodes 1 and -5
2-element Vector{Union{Nothing, Int64}}:
 1
 7

julia> params = ParamsBM(10, 0.1, ShiftNet(net.node[nodes_shifts], [3.0, -3.0],  net))
ParamsBM:
Parameters of a BM with fixed root:
mu: 10.0
Sigma2: 0.1

There are 2 shifts on the network:
──────────────────────────
  Edge Number  Shift Value
──────────────────────────
          8.0         -3.0
          1.0          3.0
──────────────────────────

julia> using Random; Random.seed!(2468); # sets the seed for reproducibility

julia> sim = rand(net, params); # simulate a dataset with shifts

julia> using DataFrames # to handle data frames

julia> dat = DataFrame(trait = sim[:tips], tipnames = sim.M.tipnames);

julia> dat = DataFrame(trait = [13.391976856737717, 9.55741491696386, 7.17703734817448, 7.889062527849697],
        tipnames = ["A","B","C","D"]) # hard-coded, to be independent of random number generator
4×2 DataFrame
 Row │ trait     tipnames 
     │ Float64   String   
─────┼────────────────────
   1 │ 13.392    A
   2 │  9.55741  B
   3 │  7.17704  C
   4 │  7.88906  D

julia> dfr_shift = descendencedataframe(net.node[nodes_shifts], net) # the regressors matching the shifts.
4×3 DataFrame
 Row │ shift_1  shift_8  tipnames 
     │ Float64  Float64  String   
─────┼────────────────────────────
   1 │     1.0      0.0  A
   2 │     0.0      0.0  B
   3 │     0.0      1.0  C
   4 │     0.0      0.6  D

julia> dfr = innerjoin(dat, dfr_shift, on=:tipnames); # join data and regressors in a single dataframe

julia> using StatsModels # for statistical model formulas

julia> fitBM = phylolm(@formula(trait ~ shift_1 + shift_8), dfr, net; reml=false) # actual fit
PhyloNetworkLinearModel

Formula: trait ~ 1 + shift_1 + shift_8

Model: Brownian motion

Parameter Estimates, using ML:
phylogenetic variance rate: 0.0112618

Coefficients:
────────────────────────────────────────────────────────────────────────
                Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
────────────────────────────────────────────────────────────────────────
(Intercept)   9.48238    0.327089  28.99    0.0220    5.32632   13.6384
shift_1       3.9096     0.46862    8.34    0.0759   -2.04479    9.86399
shift_8      -2.4179     0.422825  -5.72    0.1102   -7.7904     2.95461
────────────────────────────────────────────────────────────────────────
Log Likelihood: 1.8937302027
AIC: 4.2125395947
```

Next we illustrate the model with heterosis, aka transgressive evolution: with
a shift in the trait after successful hybridization.
First how to simulated according to this model:

```jldoctest descendence
julia> nodes_hybrids = indexin([5], [n.number for n in net.node]); # Put a shift on edges below hybrids

julia> params = ParamsBM(10, 0.1, ShiftNet(net.node[nodes_hybrids], [3.0],  net));

julia> using Random; Random.seed!(2468); # sets the seed for reproducibility

julia> sim = rand(net, params); # simulate a dataset with shifts

julia> dat = DataFrame(trait = sim[:tips], tipnames = sim.M.tipnames);
```

and next how to analyze data under a transgressive evolution model.
Below we hard-code data values for more reproducibility.

```jldoctest descendence
julia> dat = DataFrame(trait = [10.391976856737717, 9.55741491696386, 10.17703734817448, 12.689062527849698],
          tipnames = ["A","B","C","D"])
4×2 DataFrame
 Row │ trait     tipnames 
     │ Float64   String   
─────┼────────────────────
   1 │ 10.392    A
   2 │  9.55741  B
   3 │ 10.177    C
   4 │ 12.6891   D

julia> dfr_hybrid = descendencedataframe(net) # the regressors matching the hybrids (all of them)
4×3 DataFrame
 Row │ shift_6  tipnames  sum     
     │ Float64  String    Float64 
─────┼────────────────────────────
   1 │     0.0  A             0.0
   2 │     0.0  B             0.0
   3 │     0.0  C             0.0
   4 │     1.0  D             1.0

julia> dfr = innerjoin(dat, dfr_hybrid, on=:tipnames); # join data and regressors in a single dataframe

julia> using StatsModels

julia> fitBM = phylolm(@formula(trait ~ shift_6), dfr, net; reml=false) # actual fit
PhyloNetworkLinearModel

Formula: trait ~ 1 + shift_6

Model: Brownian motion

Parameter Estimates, using ML:
phylogenetic variance rate: 0.041206

Coefficients:
────────────────────────────────────────────────────────────────────────
                Coef.  Std. Error      t  Pr(>|t|)  Lower 95%  Upper 95%
────────────────────────────────────────────────────────────────────────
(Intercept)  10.064      0.277959  36.21    0.0008    8.86805   11.26
shift_6       2.72526    0.315456   8.64    0.0131    1.36796    4.08256
────────────────────────────────────────────────────────────────────────
Log Likelihood: -0.7006021946
AIC: 7.4012043891
```

# See also
[`phylolm`](@ref), [`PhyloNetworks.descendencematrix`](@extref).
"""
function descendencedataframe(
    node::Vector{Node},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    T = PN.descendencematrix(net; checkpreorder=checkpreorder)
    descendencedataframe(node, net, T)
end

function descendencedataframe(
    node::Vector{Node},
    net::HybridNetwork,
    T::MatrixTopologicalOrder
)
    ## Get the descendence matrix for tips
    T_t = T[:tips]
    ## Get the indices of the columns to keep
    ind = zeros(Int, length(node))
    for (i,nod) in enumerate(node)
        !nod.hybrid || error("Shifts on hybrid edges are not allowed")
        ii = findfirst(n -> n===nod, net.vec_node)
        isnothing(ii) && error("node number $(nod.number) (i=$i) not in preorder list")
        ind[i] = ii
    end
    ## get column names
    eNum = [getMajorParentEdgeNumber(n) for n in net.vec_node[ind]]
    function tmp_fun(x::Int)
        return(Symbol("shift_$(x)"))
    end
    df = DataFrame(T_t[:, ind], [tmp_fun(num) for num in eNum])
    df[!,:tipnames]=T.tipnames
    return(df)
end

function descendencedataframe(
    edge::Vector{Edge},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    childs = [getchild(ee) for ee in edge]
    return(descendencedataframe(childs, net; checkpreorder=checkpreorder))
end

descendencedataframe(edge::Edge, net::HybridNetwork; checkpreorder::Bool=true) = descendencedataframe([edge], net; checkpreorder=checkpreorder)
descendencedataframe(node::Node, net::HybridNetwork; checkpreorder::Bool=true) = descendencedataframe([node], net; checkpreorder=checkpreorder)

function descendencedataframe(
    net::HybridNetwork;
    checkpreorder::Bool=true,
    which=:allhybrids,
)
    if which == :allhybrids
        childs = [getchild(nn) for nn in net.hybrid] # checks that each hybrid node has a single child
        dfr = descendencedataframe(childs, net; checkpreorder=checkpreorder)
        dfr[!,:sum] = sum.(eachrow(select(dfr, Not(:tipnames), copycols=false)))
        return(dfr)
    end
    throw(ArgumentError("""`which` must be equal to `:allhybrids` for now.
    Otherwise call descendencedataframe with a vector of edges or nodes."""))
end
