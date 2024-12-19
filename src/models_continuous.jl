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

# function Base.:(==)(sh1::ShiftNet, sh2::ShiftNet)
#     isEqual(sh1.net, sh2.net) || return(false)
#     sh1.shift == sh2.shift || return(false)
#     return(true)
# end


abstract type ParamsProcess end

"""
    ParamsBM <: ParamsProcess

Type for a BM process on a network. Fields are `mu` (expectation),
`sigma2` (variance), `randomRoot` (whether the root is random, default to `false`),
and `varRoot` (if the root is random, the variance of the root, default to `NaN`).

"""
mutable struct ParamsBM <: ParamsProcess
    mu::Float64 # Ancestral value or mean
    sigma2::Float64 # variance
    randomRoot::Bool # Root is random ? default false
    varRoot::Float64 # root variance. Default NaN
    shift::Union{ShiftNet, Missing} # shifts

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

