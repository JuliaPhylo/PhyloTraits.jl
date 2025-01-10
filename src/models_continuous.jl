# types for shifts. used later in types of simulation parameters
"""
    ShiftNet

Shifts mapped to tree nodes and their (unique) parent edge on a
[`PhyloNetworks.HybridNetwork`](@extref) sorted in topological order.
Its `shift` field is a vector of shift values, one for each node,
corresponding to the shift on the parent edge of the node
(which makes sense for tree nodes only: they have a single parent edge).

Two `ShiftNet` objects on the same network can be concatened with `*`.

`ShiftNet(node::Vector{Node}, value::AbstractVector, net::HybridNetwork; checkpreorder::Bool=true)`

Constructor from a vector of nodes and associated values. The shifts are located
on the edges above the nodes provided. Warning, shifts on hybrid edges are not
allowed.

`ShiftNet(edge::Vector{Edge}, value::AbstractVector, net::HybridNetwork; checkpreorder::Bool=true)`

Constructor from a vector of edges and associated values.
Warning, shifts on hybrid edges are not allowed.

Extractors: [`getshiftedgenumber`](@ref), [`getshiftvalue`](@ref)
"""
struct ShiftNet
    shift::Matrix{Float64}
    net::HybridNetwork
end

ShiftNet(net::HybridNetwork, dim::Int) = ShiftNet(zeros(length(net.node), dim), net)
ShiftNet(net::HybridNetwork) = ShiftNet(net, 1)

# from vector of nodes
function ShiftNet(
    node::Vector{Node},
    value::AbstractMatrix,
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    n_nodes, dim = size(value)
    if length(node) != n_nodes
        error("The vector of nodes/edges and of values must have the same number or rows.")
    end
    if checkpreorder
        preorder!(net)
    end
    obj = ShiftNet(net, dim)
    for i in 1:length(node)
        !node[i].hybrid || error("Shifts on hybrid edges are not allowed")
        ind = findfirst(x -> x===node[i], net.vec_node)
        obj.shift[ind, :] .= @view value[i, :]
    end
    return(obj)
end
function ShiftNet(
    node::Vector{Node},
    value::AbstractVector,
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    return ShiftNet(node, reshape(value, (length(value), 1)), net,
                    checkpreorder = checkpreorder)
end
# from vector of edges
function ShiftNet(
    edge::Vector{Edge},
    value::Union{AbstractVector, AbstractMatrix},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    childs = [getchild(ee) for ee in edge]
    return(ShiftNet(childs, value, net; checkpreorder=checkpreorder))
end
# from a single node or single edge
ShiftNet(edge::Edge, value::Float64, net::HybridNetwork; checkpreorder::Bool=true) =
    ShiftNet([edge], [value], net; checkpreorder=checkpreorder)
ShiftNet(node::Node, value::Float64, net::HybridNetwork; checkpreorder::Bool=true) =
    ShiftNet([node], [value], net; checkpreorder=checkpreorder)
function ShiftNet(
    edge::Edge,
    value::AbstractVector{Float64},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    return ShiftNet([edge], reshape(value, (1, length(value))), net,
                    checkpreorder = checkpreorder)
end
function ShiftNet(
    node::Node,
    value::AbstractVector{Float64},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    return ShiftNet([node], reshape(value, (1, length(value))), net,
                    checkpreorder = checkpreorder)
end

"""
    shiftathybrids(value::Vector{T} where T<:Real, net::HybridNetwork; checkpreorder::Bool=true)

Construct an object [`ShiftNet`](@ref) with shifts on all the edges below
hybrid nodes, with values provided. The vector of values must have the
same length as the number of hybrids in the network.

"""
function shiftathybrids(
    value::Union{Matrix{T}, Vector{T}} where T<:Real,
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    if length(net.hybrid) != size(value, 1)
        error("You must provide as many values as the number of hybrid nodes.")
    end
    childs = [getchild(nn) for nn in net.hybrid] # checks for single child
    return(ShiftNet(childs, value, net; checkpreorder=checkpreorder))
end
shiftathybrids(value::Real, net::HybridNetwork; checkpreorder::Bool=true) =
    shiftathybrids([value], net; checkpreorder=checkpreorder)

"""
    getshiftedgenumber(shift::ShiftNet)

Get the edge numbers where the shifts are located, for an object [`ShiftNet`](@ref).
If a shift is placed at the root node with no parent edge, the edge number
of a shift is set to -1 (as if missing).
"""
function getshiftedgenumber(shift::ShiftNet)
    nodInd = getshiftrowinds(shift)
    [getMajorParentEdgeNumber(n) for n in shift.net.vec_node[nodInd]]
end

function getMajorParentEdgeNumber(n::Node)
    try
        getparentedge(n).number
    catch
        -1
    end
end

function getshiftrowinds(shift::ShiftNet)
    n, p = size(shift.shift)
    inds = zeros(Int, n)
    counter = 0
    for i = 1:n
        use_row = !all(iszero, @view shift.shift[i, :])
        if use_row
            counter += 1
            inds[counter] = i
        end
    end
    return inds[1:counter]
end
"""
    getshiftvalue(shift::ShiftNet)

Get the values of the shifts, for an object [`ShiftNet`](@ref).
"""
function getshiftvalue(shift::ShiftNet)
    rowInds = getshiftrowinds(shift)
    shift.shift[rowInds, :]
end

function shift_coeftable(shift::ShiftNet)
    sv = getshiftvalue(shift)
    if size(sv, 2) == 1
        shift_labels = ["Shift Value"]
    else
        shift_labels = ["Shift Value $i" for i = 1:size(sv, 2)]
    end
    CoefTable(hcat(getshiftedgenumber(shift), sv),
              ["Edge Number"; shift_labels],
              fill("", size(sv, 1)))
end

function Base.show(io::IO, obj::ShiftNet)
    println(io, "$(typeof(obj)):\n",
            shift_coeftable(obj))
end

function Base.:*(sh1::ShiftNet, sh2::ShiftNet)
    PN.isEqual(sh1.net, sh2.net) ||
        error("Shifts to be concatenated must be defined on the same network.")
    size(sh1.shift) == size(sh2.shift) ||
        error("Shifts to be concatenated must have the same dimensions.")
    shiftNew = zeros(size(sh1.shift))
    for i in 1:length(sh1.shift)
        if iszero(sh1.shift[i])
            shiftNew[i] = sh2.shift[i]
        elseif iszero(sh2.shift[i])
            shiftNew[i] = sh1.shift[i]
        elseif sh1.shift[i] == sh2.shift[i]
            shiftNew[i] = sh1.shift[i]
        else
            error("The two shifts matrices you provided affect the same " *
                  "trait for the same edge, so I cannot choose which one you want.")
        end
    end
    return(ShiftNet(shiftNew, sh1.net))
end

# types used for simulation, and for ancestralreconstruction

abstract type ParamsProcess end

"""
    ParamsBM <: ParamsProcess

Type for a BM process on a network. Fields are `mu` (expectation),
`sigma2` (variance), `randomRoot` (whether the root is random, default to `false`),
and `varRoot` (if the root is random, the variance of the root, default to `NaN`).
"""
mutable struct ParamsBM <: ParamsProcess
    "ancestral value (or mean)"
    mu::Float64
    "evolutionary variance"
    sigma2::Float64
    "is the root trait a random variable? false by default in many constructors"
    randomRoot::Bool
    "variance of the root trait. NaN by default in many constructors"
    varRoot::Float64
    "shifts, if any: ShiftNet object, or missing"
    shift::Union{ShiftNet, Missing}

    function ParamsBM(
        mu::Real,
        sigma2::Real,
        randomRoot::Bool,
        varRoot::Real,
        shift::Union{ShiftNet, Missing}
    )
        if !ismissing(shift) && size(shift.shift, 2) != 1
            error("ShiftNet must have only a single shift dimension.")
        end
        return new(mu, sigma2, randomRoot, varRoot, shift)
    end
end
# Constructor
ParamsBM(mu::Real, sigma2::Real) = ParamsBM(mu, sigma2, false, NaN, missing) # default values
ParamsBM(mu::Real, sigma2::Real, net::HybridNetwork) = ParamsBM(mu, sigma2, false, NaN, ShiftNet(net)) # default values
ParamsBM(mu::Real, sigma2::Real, shift::ShiftNet) = ParamsBM(mu, sigma2, false, NaN, shift) # default values

function anyShift(params::ParamsProcess)
    if ismissing(params.shift) return(false) end
    for v in params.shift.shift
        if v != 0 return(true) end
    end
    return(false)
end

process_dim(::ParamsBM) = 1

function Base.show(io::IO, obj::ParamsBM)
    disp =  "$(typeof(obj)):\n"
    pt = paramstable(obj)
    if obj.randomRoot
        disp = disp * "Parameters of a BM with random root:\n" * pt
    else
        disp = disp * "Parameters of a BM with fixed root:\n" * pt
    end
    println(io, disp)
end

function paramstable(obj::ParamsBM)
    disp = "mu: $(obj.mu)\nSigma2: $(obj.sigma2)"
    if obj.randomRoot
        disp = disp * "\nvarRoot: $(obj.varRoot)"
    end
    if anyShift(obj)
        disp = disp * "\n\nThere are $(length(getshiftvalue(obj.shift))) shifts on the network:\n"
        disp = disp * "$(shift_coeftable(obj.shift))"
    end
    return(disp)
end

"""
    ParamsMultiBM <: ParamsProcess

Type for a multivariate Brownian diffusion (MBD) process on a network. Fields are `mu` (expectation),
`sigma` (covariance matrix), `randomRoot` (whether the root is random, default to `false`),
`varRoot` (if the root is random, the covariance matrix of the root, default to `[NaN]`),
`shift` (a ShiftNet type, default to `missing`),
and `L` (the lower triangular of the cholesky decomposition of `sigma`, computed automatically)

# examples and constructors

```jldoctest
julia> ParamsMultiBM([1.0, -0.5], [2.0 0.3; 0.3 1.0]) # no shifts
ParamsMultiBM:
Parameters of a MBD with fixed root:
mu: [1.0, -0.5]
Sigma: [2.0 0.3; 0.3 1.0]

julia> net = readnewick("((A:1,B:1):1,C:2);");

julia> shifts = ShiftNet(net.node[2], [-1.0, 2.0], net);

julia> ParamsMultiBM([1.0, -0.5], [2.0 0.3; 0.3 1.0], shifts) # with shifts
ParamsMultiBM:
Parameters of a MBD with fixed root:
mu: [1.0, -0.5]
Sigma: [2.0 0.3; 0.3 1.0]

There are 2 shifts on the network:
───────────────────────────────────────────
  Edge Number  Shift Value 1  Shift Value 2
───────────────────────────────────────────
          2.0           -1.0            2.0
───────────────────────────────────────────

```
"""
mutable struct ParamsMultiBM <: ParamsProcess
    mu::Vector{Float64}
    sigma::Matrix{Float64}
    randomRoot::Bool
    varRoot::Matrix{Float64}
    shift::Union{ShiftNet, Missing}
    L::LowerTriangular{Float64}

    function ParamsMultiBM(
        mu::AbstractArray{Float64, 1},
        sigma::AbstractArray{Float64, 2},
        randomRoot::Bool,
        varRoot::AbstractArray{Float64, 2},
        shift::Union{ShiftNet, Missing},
        L::LowerTriangular{Float64}
    )
        dim = length(mu)
        if size(sigma) != (dim, dim)
            error("The mean and variance do must have conforming dimensions.")
        end
        if randomRoot && size(sigma) != size(varRoot)
            error("The root variance and process variance must have the same dimensions.")
        end
        if !ismissing(shift) && size(shift.shift, 2) != dim
            error("The ShiftNet and diffusion process must have the same dimensions.")
        end
        return new(mu, sigma, randomRoot, varRoot, shift, L)
    end
end

ParamsMultiBM(mu::AbstractArray{Float64, 1}, sigma::AbstractArray{Float64, 2}) =
    ParamsMultiBM(mu, sigma, false, Diagonal([NaN]), missing, cholesky(sigma).L)
function ParamsMultiBM(
    mu::AbstractArray{Float64, 1},
    sigma::AbstractArray{Float64, 2},
    shift::ShiftNet
)
    return ParamsMultiBM(mu, sigma, false, Diagonal([NaN]), shift, cholesky(sigma).L)
end
function ParamsMultiBM(
    mu::AbstractArray{Float64, 1},
    sigma::AbstractArray{Float64, 2},
    net::HybridNetwork
)
    return ParamsMultiBM(mu, sigma, ShiftNet(net, length(mu)))
end

process_dim(params::ParamsMultiBM) = length(params.mu)

function Base.show(io::IO, obj::ParamsMultiBM)
    disp =  "$(typeof(obj)):\n"
    pt = paramstable(obj)
    if obj.randomRoot
        disp = disp * "Parameters of a MBD with random root:\n" * pt
    else
        disp = disp * "Parameters of a MBD with fixed root:\n" * pt
    end
    println(io, disp)
end

function paramstable(obj::ParamsMultiBM)
    disp = "mu: $(obj.mu)\nSigma: $(obj.sigma)"
    if obj.randomRoot
        disp = disp * "\nvarRoot: $(obj.varRoot)"
    end
    if anyShift(obj)
        disp = disp * "\n\nThere are $(length(getshiftvalue(obj.shift))) shifts on the network:\n"
        disp = disp * "$(shift_coeftable(obj.shift))"
    end
    return(disp)
end


"""
    WithinSpeciesCTM

Type to fit models accounting for within-species variation, including
measurement error, genetic variation between individuals, plasticity,
environmental variation etc.
CTM stands for "continuous trait model". Contains the estimated variance components
(between-species phylogenetic variance rate and within-species variance)
and output from the `NLopt` optimization used in the estimation.

## fields

- `wsp_var`: intra/within-species variance.
- `bsp_var`: inter/between-species variance-rate.
- `wsp_ninv`: vector of the inverse sample-sizes (e.g. [1/n₁, ..., 1/nₖ], where
  data from k species was used to fit the model and nᵢ is the no. of observations
  for the ith species).
- `rss`: within-species sum of squares
- `optsum`: an [`OptSummary`](@ref) object.
"""
struct WithinSpeciesCTM
    "within-species variance η*σ², assumes Normal distribution"
    wsp_var::Vector{Float64} # vector to make it mutable
    "between-species variance rate σ², such as from Brownian motion"
    bsp_var::Vector{Float64}
    "inverse sample sizes (or precision): 1/(no. of individuals) within each species"
    wsp_ninv::Vector{Float64}
    "within-species sum of squares"
    rss::Float64
    "NLopt & NLopt summary object"
    optsum::OptSummary
end

"""
    ContinuousTraitEM

Abstract type for evolutionary models for continuous traits, using a continuous-time
stochastic process on a phylogeny.

For subtypes, see [`BM`](@ref), [`PagelLambda`](@ref), [`ScalingHybrid`](@ref).

Each of these subtypes/models has the field `lambda`, whose default value is 1.0.
However, the interpretation of this field differs across models.
"""
abstract type ContinuousTraitEM end

# current concrete subtypes: BM, PagelLambda, ScalingHybrid
# possible future additions: OU (Ornstein-Uhlenbeck)?
"""
    BM(λ)

Brownian Motion, subtype of [`ContinuousTraitEM`](@ref), to model the population mean
of a trait (or of the residuals from a linear model). Under the BM model,
the population (or species) means have a multivariate normal distribution with
covariance matrix = σ²λV, where σ² is the between-species
variance-rate (to be estimated), and the matrix V is obtained from
[`PhyloNetworks.sharedpathmatrix`](@extref)`(net)[:tips]`.

λ is set to 1 by default, and is immutable.
In future versions, λ may be used to control the scale for σ².

On a tree, V is the length of shared ancestry.
On a network, the BM model assumes that the trait at a hybrid node
is the weighted average of its immediate parents (plus possibly a fixed shift).
The weights are the proportion of genes inherited from each parent:
the γ parameters of hybrid edges.
"""
struct BM <: ContinuousTraitEM
    lambda::Float64 # immutable
end
BM() = BM(1.0)
evomodelname(::BM) = "Brownian motion"

"""
    PagelLambda(λ)

Pagel's λ model, subtype of [`ContinuousTraitEM`](@ref), with covariance matrix σ²V(λ).
σ² is the between-species variance-rate (to be estimated), and V(λ) = λV + (1-λ)T,
where V is the covariance under a Brownian motion [`BM`](@ref) and T is a diagonal
matrix containing the total branch length elapsed from the root to each leaf (if
the phylogeny is a tree, or more generally if the network is time consistent: the
time from the root to a given node does not depend on the path).

λ ∈ [0,1] is mutable and may be optimized. It is a measure of phylogenetic
signal, that is, how important the given network is for explaining variation in
the response. When λ=1, the `PagelLambda` model reduces to the `BM` model.
"""
mutable struct PagelLambda <: ContinuousTraitEM
    lambda::Float64 # mutable: can be optimized
end
PagelLambda() = PagelLambda(1.0)
evomodelname(::PagelLambda) = "Pagel's lambda"

"""
    ScalingHybrid(λ)

Scaling Hybrid model, subtype of [`ContinuousTraitEM`](@ref), with covariance matrix
σ²V(N(λ)). σ² is the between-species variance-rate (to be estimated),
V(N) is the Brownian motion [`BM`](@ref) covariance obtained from network N,
and N(λ) is a obtained from the input network by rescaling the inheritance
parameter γ of all minor edges by the same λ: a minor edge has its original γ
changed to λγ, using the same λ at all reticulations.
Note that for a major edge with original inheritance γ, the partner minor edge
has inheritance γ_minor = 1-γ, so the major edge's inheritance is changed to
1-λγ_minor = λγ+1-λ.

For more information: see Bastide (2017) dissertation, section 4.3.2 p.175,
available at https://tel.archives-ouvertes.fr/tel-01629648

λ ∈ [0,1] is mutable and may be optimized. It is a measure of how important the
reticulations are for explaining variation in the response.
When λ=1, the `ScalingHybrid` model reduces to the `BM` model.
"""
mutable struct ScalingHybrid <: ContinuousTraitEM
    lambda::Float64
end
ScalingHybrid() = ScalingHybrid(1.0)
evomodelname(::ScalingHybrid) = "Lambda's scaling hybrid"

###############################################################################
##     phylogenetic network regression
###############################################################################

"""
    PhyloNetworkLinearModel <: GLM.LinPredModel

Phylogenetic linear model representation.

## Fields

`lm`, `V`, `Vy`, `RL`, `Y`, `X`, `logdetVy`, `reml`, `ind`, `nonmissing`,
`evomodel`, `model_within` and `formula`.
The following syntax pattern can be used to get more information on a specific field:
e.g. to find out about the `lm` field, do `?PhyloNetworkLinearModel.lm`.

## Methods applied to fitted models

The following StatsAPI / StatsBase functions can be applied:
`coef`, `nobs`, `vcov`, `stderror`, `confint`, `coeftable`, `dof_residual`, `dof`, `deviance`,
`residuals`, `response`, `predict`, `loglikelihood`, `nulldeviance`, `nullloglikelihood`,
`r2`, `adjr2`, `aic`, `aicc`, `bic`, `ftest`, `lrtest` etc.

The estimated variance-rate and estimated mean of the species-level trait model
(see [`ContinuousTraitEM`](@ref)) can be retrieved using [`sigma2_phylo`](@ref)
and [`mu_phylo`](@ref) respectively.

If relevant, the estimated individual-level/within-species variance can be retrieved
using [`sigma2_within`](@ref).

The optimized λ parameter for Pagel's λ model (see [`PagelLambda`](@ref)) can
be retrieved using [`lambda_estim`](@ref).

An ancestral state reconstruction can be performed using [`ancestralreconstruction`](@ref).

## Within-species variation

The true species/population means for the response trait/variable (or the residuals:
conditional on the predictors) are jointly modeled as 𝒩(·, σ²ₛV) where V depends on
the trait model (see [`ContinuousTraitEM`](@ref)) and on the species network.
σ²ₛ is the between-species variance-rate.

Within-species variation is modeled by assuming that the individual-level
responses are iid 𝒩(0, σ²ₑ) about the true species means, so that the
species-level sample means (conditional on the predictors) are jointly modeled
as 𝒩(·, σ²ₛV + σ²ₑD⁻¹), where σ²ₑ is the within-species variance and D⁻¹ is a
diagonal matrix whose entries are the inverse sample-sizes (see [`WithinSpeciesCTM`](@ref)).

Although the above two models can be expressed in terms of a joint distribution
for the species-level sample means (or residuals conditional on the predictors),
more data are required to fit a model accounting for within-species variation,
that is, a model recognizing that the sample means are estimates of the true
population means. To fit a model *without* within-species variation, data on the
species means are sufficient. To fit a model *with* within-species variation,
we need to have the species means and the standard deviations of the response
variable for each species.

`phylolm` can fit a model with within-species variation either from
species-level statistics ("mean response" and "standard deviation in response")
or from individual-level data (in which case the species-level statistics are
computed internally). See [`phylolm`](@ref) for more details on these two
input choices.

In the object, `obj.Y` and `obj.X` are the observed species means.
`predict`, `residuals` and `response` return the values at the species level.
"""
mutable struct PhyloNetworkLinearModel <: GLM.LinPredModel
    "lm: a GLM.LinearModel object, fitted on the cholesky-tranformed problem"
    lm::GLM.LinearModel # result of a lm on a matrix
    "V: a MatrixTopologicalOrder object of the network-induced correlations"
    V::MatrixTopologicalOrder
    "Vy: the sub matrix corresponding to the tips and actually used for the correction"
    Vy::Matrix
    """RL: a LowerTriangular matrix, the lower Cholesky factor of Vy=RL*RL'
    obtained with `cholesky(Vy).L`. The data stored in `lm` are RL⁻¹Y and RL⁻¹X.
    """
    RL::LowerTriangular
    "Y: the vector of data"
    Y::Vector
    "X: the matrix of regressors"
    X::Matrix
    "logdetVy: the log-determinant of Vy"
    logdetVy::Float64
    "criterion: REML if reml is true, ML otherwise"
    reml::Bool
    "ind: vector matching the tips of the network against the names of the dataframe provided. 0 if the match could not be performed."
    ind::Vector{Int}
    "nonmissing: vector indicating which tips have non-missing data"
    nonmissing::BitArray{1}
    "evomodel: the model used for the fit"
    evomodel::ContinuousTraitEM
    # ContinuousTraitEM is abstract: not efficient. parametrize PhyloNetworkLinearModel?
    # but the types for Vy, Y and X are also abstract.
    "model_within: the model used for within-species variation (if needed)"
    model_within::Union{Nothing, WithinSpeciesCTM}
    "formula: a StatsModels.FormulaTerm formula"
    formula::Union{StatsModels.FormulaTerm,Nothing}
end

# default model_within=nothing
PhyloNetworkLinearModel(lm,  V,Vy,RL,Y,X,logdetVy, reml,ind,nonmissing, model) =
  PhyloNetworkLinearModel(lm,V,Vy,RL,Y,X,logdetVy, reml,ind,nonmissing, model,nothing,nothing)
# default formula=nothing
PhyloNetworkLinearModel(lm,V,Vy,RL,Y,X,logdetVy,reml,ind,nonmissing,model,model_within) =
PhyloNetworkLinearModel(lm,V,Vy,RL,Y,X,logdetVy, reml,ind,nonmissing,model,model_within,nothing)

#= Methods for PhyloNetworkLinearModel objects
un-changed quantities:
- regression coefficients
- formula
- hasintercept
As PhyloNetworkLinearModel <: GLM.LinPredModel, these are automatically defined:
aic, aicc, bic (so long as there is a valid method for loglikelihood, nobs etc.)
=#
StatsAPI.coef(m::PhyloNetworkLinearModel) = coef(m.lm)
StatsAPI.coefnames(m::PhyloNetworkLinearModel) =
    (m.formula === nothing ? ["x$i" for i in 1:length(coef(m))] : coefnames(formula(m).rhs))
StatsModels.formula(obj::PhyloNetworkLinearModel) = obj.formula
StatsModels.hasintercept(m::PhyloNetworkLinearModel) =
    any(i -> all(==(1), view(m.X , :, i)), 1:size(m.X, 2))

"""
    StatsBase.deviance(m::PhyloNetworkLinearModel)

-2 loglikelihood of the fitted model. See also  [`loglikelihood`](@ref).

Note: this is not the residual-sum-of-squares deviance as output by GLM,
such as one would get with `deviance(m.model)`.
"""
function StatsBase.deviance(m::PhyloNetworkLinearModel, ::Val{false}=Val(false))
    -2*loglikelihood(m)
end
"""
    StatsBase.deviance(m::PhyloNetworkLinearModel, Val(true))

Residual sum of squares with metric V, the estimated phylogenetic covariance,
if the model is appropriate.
"""
function StatsBase.deviance(m::PhyloNetworkLinearModel, ::Val{true})
    isnothing(m.model_within) ||
        error("deviance measured as SSR not implemented for within-species variation")
    deviance(m.lm)
end
"""
    StatsBase.nobs(m::PhyloNetworkLinearModel)

Number of observations: number of species with data, if the model assumes
known species means, and number of individuals with data, if the model
accounts for within-species variation.
"""
function StatsBase.nobs(m::PhyloNetworkLinearModel)
    if isnothing(m.model_within)
        return nobs(m.lm)
    else
        return sum(1.0 ./ m.model_within.wsp_ninv)
    end
end

# degrees of freedom for residuals: at the species level, for coefficients of
# the phylogenetic regression, assuming known co-variance between species means.
# used for F and T degrees of freedom, instead of more conservative Z
StatsAPI.dof_residual(m::PhyloNetworkLinearModel) =  nobs(m.lm) - length(coef(m))

# degrees of freedom consumed by the species-level model
function StatsAPI.dof(m::PhyloNetworkLinearModel)
    res = length(coef(m)) + 1 # +1: phylogenetic variance
    if any(typeof(m.evomodel) .== [PagelLambda, ScalingHybrid])
        res += 1 # lambda is one parameter
    end
    if !isnothing(m.model_within)
        res += 1 # within-species variance
    end
    return res
end

"""
    vcov(m::PhyloNetworkLinearModel)

Return the variance-covariance matrix of the coefficient estimates.

For the continuous trait evolutionary models currently implemented, species-level
mean response (conditional on the predictors), Y|X is modeled as:

1. Y|X ∼ 𝒩(Xβ, σ²ₛV) for models assuming known species mean values (no within-species variation)
2. Y|X ∼ 𝒩(Xβ, σ²ₛV + σ²ₑD⁻¹) for models with information from multiple individuals
   and assuming within-species variation

The matrix V is inferred from the phylogeny, but may also depend on additional
parameters to be estimated (e.g. `lambda` for Pagel's Lambda model). See
[`ContinuousTraitEM`](@ref), [`PhyloNetworkLinearModel`](@ref) for more details.

If (1), then return σ²ₛ(X'V⁻¹X)⁻¹, where σ²ₛ is estimated with REML, even if
the model was fitted with `reml=false`.
This follows the conventions of [`nlme::gls`](https://www.rdocumentation.org/packages/nlme/versions/3.1-152)
and [`stats::glm`](https://www.rdocumentation.org/packages/stats/versions/3.6.2) in R.

If (2), then return σ²ₛ(X'W⁻¹X)⁻¹, where W = V+(σ²ₑ/σ²ₛ)D⁻¹ is estimated, and
σ²ₛ & σₑ are the estimates obtained with ML or REML, depending on the `reml`
option used to fit the model `m`. This follows the convention
of [`MixedModels.fit`](https://juliastats.org/MixedModels.jl/stable/) in Julia.
"""
function StatsBase.vcov(m::PhyloNetworkLinearModel)
    # GLM.vcov(x::LinPredModel) = rmul!(invchol(x.pp), dispersion(x, true))
    # GLM.dispersion (sqrt=true): sqrt(sum of working residuals / dof_residual): forces "REML"
    (isnothing(m.model_within) ? vcov(m.lm) :
                                 rmul!(GLM.invchol(m.lm.pp), sigma2_phylo(m)) )
end

"""
    stderror(m::PhyloNetworkLinearModel)

Return the standard errors of the coefficient estimates. See [`vcov`](@ref)
for related information on how these are computed.
"""
StatsBase.stderror(m::PhyloNetworkLinearModel) = sqrt.(diag(vcov(m)))
# confidence intervals for coefficients: GLM uses normal quantiles
# Based on: https://github.com/JuliaStats/GLM.jl/blob/d1ccc9abcc9c7ca6f640c13ff535ee8383e8f808/src/lm.jl#L240-L243

"""
    confint(m::PhyloNetworkLinearModel; level::Real=0.95)

Return confidence intervals for coefficients, with confidence level `level`,
based on the t-distribution whose degree of freedom is determined by the
number of species (as returned by `dof_residual`)
"""
function StatsBase.confint(m::PhyloNetworkLinearModel; level::Real=0.95)
    hcat(coef(m),coef(m)) + stderror(m) *
    quantile(TDist(dof_residual(m)), (1. - level)/2.) * [1. -1.]
end
# Table of estimated coefficients, standard errors, t-values, p-values, CIs
# Based on: https://github.com/JuliaStats/GLM.jl/blob/d1ccc9abcc9c7ca6f640c13ff535ee8383e8f808/src/lm.jl#L193-L203
"""
    coeftable(m::PhyloNetworkLinearModel; level::Real=0.95)

Return coefficient estimates, standard errors, t-values, p-values, and t-intervals
as a `StatsBase.CoefTable`.
"""
function StatsAPI.coeftable(m::PhyloNetworkLinearModel; level::Real=0.95)
    n_coef = size(m.lm.pp.X, 2) # no. of predictors
    if n_coef == 0
        return CoefTable([0], ["Fixed Value"], ["(Intercept)"])
    end
    cc = coef(m)
    se = stderror(m)
    tt = cc ./ se
    p = ccdf.(Ref(FDist(1, dof_residual(m))), abs2.(tt))
    ci = se*quantile(TDist(dof_residual(m)), (1-level)/2)
    levstr = isinteger(level*100) ? string(Integer(level*100)) : string(level*100)
    cn = StatsModels.vectorize(coefnames(m))
    CoefTable(hcat(cc,se,tt,p,cc+ci,cc-ci),
              ["Coef.","Std. Error","t","Pr(>|t|)","Lower $levstr%","Upper $levstr%"],
              cn, 4, 3)
end
#= changed quantities: those using residuals rescaled by cholesky of variance
   between tips
=#
StatsAPI.residuals(m::PhyloNetworkLinearModel) = m.RL * residuals(m.lm)
# Tip data: m.Y is different from response(m.lm)
#       and m.X is different from modelmatrix(m.lm)
StatsAPI.response(m::PhyloNetworkLinearModel) = m.Y
StatsAPI.modelmatrix(m::PhyloNetworkLinearModel) = m.X
# Predicted values at the tips
# (rescaled by cholesky of tips variances)
StatsAPI.predict(m::PhyloNetworkLinearModel) = m.RL * predict(m.lm)

"""
    loglikelihood(m::PhyloNetworkLinearModel)

Log likelihood or log restricted likelihood (REML) depending on `m.reml`,
of the fitted model.

For models with no within-species variation, the likelihood (or REML) is
calculated based on the joint density for species-level mean responses.

For within-species variation models, the likelihood is calculated based on the joint
density for individual-level responses. This can be calculated from individual-level
data, but also by providing species-level means and standard deviations which is
accepted by [`phylolm`](@ref).

**Warning**: many summaries are based on the species-level model, like
"dof_residual", "residuals", "predict" or "deviance".
So `deviance` is innapropriate to compare models with within-species variation.
Use `loglikelihood` to compare models based on data at the individual level.

Reminder: do not compare ML or REML values across models fit on different data.
Do not compare REML values across models that do not have the same predictors
(fixed effects): use ML instead, for that purpose.
"""
function StatsAPI.loglikelihood(m::PhyloNetworkLinearModel)
    linmod = m.lm
    if isnothing(m.model_within) # not a msrerr model
        n = (m.reml ? dof_residual(linmod) : nobs(linmod) )
        σ² = deviance(linmod)/n
        ll =  - n * (1. + log2π + log(σ²))/2 - m.logdetVy/2
    else # if msrerr model, return loglikelihood of individual-level data
        modwsp = m.model_within
        ntot = sum(1.0 ./ modwsp.wsp_ninv) # total number of individuals
        nsp = nobs(linmod)                 # number of species
        ncoef = length(coef(linmod))
        bdof = (m.reml ? nsp - ncoef : nsp )
        wdof = ntot - nsp
        N = wdof + bdof # ntot or ntot - ncoef
        σ²  = modwsp.bsp_var[1]
        σw² = modwsp.wsp_var[1]
        ll = sum(log.(modwsp.wsp_ninv)) -
             (N + N * log2π + bdof * log(σ²) + wdof * log(σw²) + m.logdetVy)
        ll /= 2
    end
    if m.reml
        ll -= logdet(linmod.pp.chol)/2 # -1/2 log|X'Vm^{-1}X|
    end
    return ll
end
# remark: this is *not* just the null deviance of the cholesky regression
"""
    StatsBase.nulldeviance(m::PhyloNetworkLinearModel)

For appropriate phylogenetic linear models, the deviance of the null model 
is the total sum of square with respect to the metric V,
the estimated phylogenetic covariance matrix.
"""
function StatsBase.nulldeviance(m::PhyloNetworkLinearModel)
    isnothing(m.model_within) ||
        error("""null loglik / deviance not implemented for within-species variation (mixed model):
        please fit the model with an intercept only instead.""")
    ro = response(m.lm) 
    if hasintercept(m)
        vo = ones(length(m.Y), 1)
        vo = m.RL \ vo
        bo = inv(vo'*vo)*vo'*ro
        ro = ro - vo*bo
    end
    return sum(ro.^2)
end
# loglikelihood of null model (intercept-only): same remark as above
function StatsBase.nullloglikelihood(m::PhyloNetworkLinearModel)
    nulldev = nulldeviance(m) # error & exit if within-species variation
    m.reml && @warn "ML null loglik: do not compare with REML on model with predictors"
    n = length(m.Y)
    return -n/2 * (log(2*pi * nulldev/n) + 1) - 1/2 * m.logdetVy
end
# coefficient of determination (1 - SS_res/SS_null)
# copied from GLM.jl/src/lm.jl, line 139
function StatsBase.r2(m::PhyloNetworkLinearModel)
    isnothing(m.model_within) ||
        error("r2 and adjusted r2 not implemented for within-species variation (mixed model)")
    1 - deviance(m, Val(true))/nulldeviance(m)
end
# adjusted coefficient of determination
# in GLM.jl/src/lm.jl: p = dof-1 and dof(x::LinearModel) = length(coef(x))+1, +1 for the dispersion parameter
function StatsBase.adjr2(obj::PhyloNetworkLinearModel)
    n = nobs(obj)
    # dof() includes the dispersion parameter sigma2, and lambda if relevant
    p = dof(obj)-1 # like in GLM
    # one could argue to use this: p = length(coef(obj)), to ignore lambda or other parameters
    1 - (1 - r2(obj))*(n-1)/(n-p)
end

# new quantities specific to *phylogenetic* linear models
"""
    sigma2_phylo(m::PhyloNetworkLinearModel)

Estimated between-species variance-rate for a fitted object.
"""
function sigma2_phylo(m::PhyloNetworkLinearModel)
    linmod = m.lm
    if isnothing(m.model_within)
        n = (m.reml ? dof_residual(linmod) : nobs(linmod) )
        σ² = deviance(linmod)/n
    else
        σ²  = m.model_within.bsp_var[1]
    end
    return σ²
end
"""
    sigma2_within(m::PhyloNetworkLinearModel)

Estimated within-species variance for a fitted object.
"""
sigma2_within(m::PhyloNetworkLinearModel) = (isnothing(m.model_within) ? nothing : m.model_within.wsp_var[1])
# ML estimate for ancestral state of the BM
"""
    mu_phylo(m::PhyloNetworkLinearModel)

Estimated root value for a fitted object.
"""
function mu_phylo(m::PhyloNetworkLinearModel)
    if m.formula === nothing 
        @warn """You fitted the data against a custom matrix, so I have no way
        to know which column is your intercept (column of ones).
        I am using the first coefficient for ancestral mean mu by convention,
        but that might not be what you are looking for."""
        size(m.lm.pp.X,2) == 0 && return 0
    elseif m.formula.rhs.terms[1] != StatsModels.InterceptTerm{true}()
        error("The fit was done without intercept, so I cannot estimate mu")
    end
    return coef(m)[1]
end

"""
    lambda(m::PhyloNetworkLinearModel)
    lambda(m::ContinuousTraitEM)

Value assigned to the lambda parameter, if appropriate.
"""
lambda(m::PhyloNetworkLinearModel) = lambda(m.evomodel)
lambda(m::Union{BM,PagelLambda,ScalingHybrid}) = m.lambda

"""
    lambda!(m::PhyloNetworkLinearModel, newlambda)
    lambda!(m::ContinuousTraitEM, newlambda)

Assign a new value to the lambda parameter.
"""
lambda!(m::PhyloNetworkLinearModel, lambda_new) = lambda!(m.evomodel, lambda_new)
lambda!(m::Union{BM,PagelLambda,ScalingHybrid}, lambda_new::Real) = (m.lambda = lambda_new)

"""
    lambda_estim(m::PhyloNetworkLinearModel)

Estimated lambda parameter for a fitted object.
"""
lambda_estim(m::PhyloNetworkLinearModel) = lambda(m)

# print & show results
function paramstable(m::PhyloNetworkLinearModel)
    Sig = sigma2_phylo(m)
    res = "phylogenetic variance rate: " * @sprintf("%.6g", Sig)
    if any(typeof(m.evomodel) .== [PagelLambda, ScalingHybrid])
        Lamb = lambda_estim(m)
        res = res*"\nLambda: " * @sprintf("%.6g", Lamb)
    end
    mw = m.model_within
    if !isnothing(mw)
        res = res*"\nwithin-species variance: " * @sprintf("%.6g", mw.wsp_var[1])
    end
    return(res)
end
# see also Base.show in
# https://github.com/JuliaStats/StatsModels.jl/blob/master/src/statsmodel.jl
function Base.show(io::IO, obj::PhyloNetworkLinearModel)
    ct = coeftable(obj)
    println(io, "$(typeof(obj))")
    if !(obj.formula === nothing)
        print(io, "\nFormula: ")
        println(io, string(obj.formula)) # formula
    end
    println(io)
    println(io, "Model: $(evomodelname(obj.evomodel))")
    println(io)
    println(io,"Parameter Estimates, using ", (obj.reml ? "REML" : "ML"),":")
    println(io, paramstable(obj))
    println(io)
    println(io,"Coefficients:")
    show(io, ct)
    println(io)
    println(io, "Log Likelihood: "*"$(round(loglikelihood(obj), digits=10))")
    println(io, "AIC: "*"$(round(aic(obj), digits=10))")
end
