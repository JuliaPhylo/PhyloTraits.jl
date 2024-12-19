
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
mu: 10
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

# Default
ShiftNet(net::HybridNetwork, dim::Int) = ShiftNet(zeros(length(net.node), dim), net)
ShiftNet(net::HybridNetwork) = ShiftNet(net, 1)

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

# Construct from edges and values
function ShiftNet(
    edge::Vector{Edge},
    value::Union{AbstractVector, AbstractMatrix},
    net::HybridNetwork;
    checkpreorder::Bool=true
)
    childs = [getchild(ee) for ee in edge]
    return(ShiftNet(childs, value, net; checkpreorder=checkpreorder))
end

ShiftNet(edge::Edge, value::Float64, net::HybridNetwork; checkpreorder::Bool=true) = ShiftNet([edge], [value], net; checkpreorder=checkpreorder)
ShiftNet(node::Node, value::Float64, net::HybridNetwork; checkpreorder::Bool=true) = ShiftNet([node], [value], net; checkpreorder=checkpreorder)

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
    shiftHybrid(value::Vector{T} where T<:Real, net::HybridNetwork; checkpreorder::Bool=true)

Construct an object [`ShiftNet`](@ref) with shifts on all the edges below
hybrid nodes, with values provided. The vector of values must have the
same length as the number of hybrids in the network.

"""
function shiftHybrid(
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
shiftHybrid(value::Real, net::HybridNetwork; checkpreorder::Bool=true) = shiftHybrid([value], net; checkpreorder=checkpreorder)

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

function shiftTable(shift::ShiftNet)
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
            shiftTable(obj))
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

###################################################
# types to hold parameters for evolutionary process
# like scalar BM, multivariate BM, OU?

abstract type ParamsProcess end

"""
    ParamsBM <: ParamsProcess

Type for a BM process on a network. Fields are `mu` (expectation),
`sigma2` (variance), `randomRoot` (whether the root is random, default to `false`),
and `varRoot` (if the root is random, the variance of the root, default to `NaN`).

"""
mutable struct ParamsBM <: ParamsProcess
    mu::Real # Ancestral value or mean
    sigma2::Real # variance
    randomRoot::Bool # Root is random ? default false
    varRoot::Real # root variance. Default NaN
    shift::Union{ShiftNet, Missing} # shifts

    function ParamsBM(mu::Real,
                      sigma2::Real,
                      randomRoot::Bool,
                      varRoot::Real,
                      shift::Union{ShiftNet, Missing})
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
        disp = disp * "$(shiftTable(obj.shift))"
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

# Constructors
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
    mu::AbstractArray{Float64, 1}
    sigma::AbstractArray{Float64, 2}
    randomRoot::Bool
    varRoot::AbstractArray{Float64, 2}
    shift::Union{ShiftNet, Missing}
    L::LowerTriangular{Float64}

    function ParamsMultiBM(mu::AbstractArray{Float64, 1},
                           sigma::AbstractArray{Float64, 2},
                           randomRoot::Bool,
                           varRoot::AbstractArray{Float64, 2},
                           shift::Union{ShiftNet, Missing},
                           L::LowerTriangular{Float64})
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

ParamsMultiBM(mu::AbstractArray{Float64, 1},
              sigma::AbstractArray{Float64, 2}) =
        ParamsMultiBM(mu, sigma, false, Diagonal([NaN]), missing, cholesky(sigma).L)

function ParamsMultiBM(mu::AbstractArray{Float64, 1},
                       sigma::AbstractArray{Float64, 2},
                       shift::ShiftNet)
    ParamsMultiBM(mu, sigma, false, Diagonal([NaN]), shift, cholesky(sigma).L)
end

function ParamsMultiBM(mu::AbstractArray{Float64, 1},
                       sigma::AbstractArray{Float64, 2},
                       net::HybridNetwork)
    ParamsMultiBM(mu, sigma, ShiftNet(net, length(mu)))
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
        disp = disp * "$(shiftTable(obj.shift))"
    end
    return(disp)
end


function partitionMBDMatrix(M::Matrix{Float64}, dim::Int)
    means = @view M[1:dim, :]
    vals = @view M[(dim + 1):(2 * dim), :]
    return means, vals
end

###############################################################################
## Type for models with within-species variation (including measurement error)

"""
    WithinSpeciesCTM

CTM stands for "continuous trait model". Contains the estimated variance components
(between-species phylogenetic variance rate and within-species variance)
and output from the `NLopt` optimization used in the estimation.

## Fields

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
function phylolm(f::StatsModels.FormulaTerm,
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
                y_mean_std::Bool=false)
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

### Methods on type phyloNetworkRegression

## Un-changed Quantities
# Coefficients of the regression
StatsAPI.coef(m::PhyloNetworkLinearModel) = coef(m.lm)
StatsAPI.coefnames(m::PhyloNetworkLinearModel) =
    (m.formula === nothing ? ["x$i" for i in 1:length(coef(m))] : coefnames(formula(m).rhs))
StatsModels.formula(obj::PhyloNetworkLinearModel) = obj.formula

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
# confidence Intervals for coefficients: GLM uses normal quantiles
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

## Changed Quantities
# Compute the residuals
# (Rescaled by cholesky of variance between tips)
StatsAPI.residuals(m::PhyloNetworkLinearModel) = m.RL * residuals(m.lm)
# Tip data: m.Y is different from response(m.lm)
#       and m.X is different from modelmatrix(m.lm)
StatsAPI.response(m::PhyloNetworkLinearModel) = m.Y
StatsAPI.modelmatrix(m::PhyloNetworkLinearModel) = m.X
# Predicted values at the tips
# (rescaled by cholesky of tips variances)
StatsAPI.predict(m::PhyloNetworkLinearModel) = m.RL * predict(m.lm)

# log likelihood of the fitted linear model
"""
    loglikelihood(m::PhyloNetworkLinearModel)

Log likelihood, or log restricted likelihood (REML), depending on `m.reml`.

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

# REMARK Not just the null deviance of the cholesky regression
# Might be something better to do than this, though.
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
StatsModels.hasintercept(m::PhyloNetworkLinearModel) = any(i -> all(==(1), view(m.X , :, i)), 1:size(m.X, 2))
# Null Log likelihood (null model with only the intercept)
# Same remark
function StatsBase.nullloglikelihood(m::PhyloNetworkLinearModel)
    nulldev = nulldeviance(m) # error & exit if within-species variation
    m.reml && @warn "ML null loglik: do not compare with REML on model with predictors"
    n = length(m.Y)
    return -n/2 * (log(2*pi * nulldev/n) + 1) - 1/2 * m.logdetVy
end
# coefficient of determination (1 - SS_res/SS_null)
# Copied from GLM.jl/src/lm.jl, line 139
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

## REMARK
# As PhyloNetworkLinearModel <: GLM.LinPredModel, the following functions are automatically defined:
# aic, aicc, bic

## New quantities
# ML estimate for variance of the BM
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

### Print the results
# Variance
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
# For DataFrameModel. see also Base.show in
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
function phylolm_wsp(X::Matrix, Y::Vector, V::MatrixTopologicalOrder,
        reml::Bool, nonmissing::BitArray{1}, ind::Vector{Int},
        ftolRel::AbstractFloat,
        xtolRel::AbstractFloat,
        ftolAbs::AbstractFloat,
        xtolAbs::AbstractFloat,
        counts::Union{Nothing, Vector},
        ySD::Union{Nothing, Vector})
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
function withinsp_varianceratio(X::Matrix, Y::Vector, V::Matrix, reml::Bool,
        ftolRel::AbstractFloat, xtolRel::AbstractFloat,
        ftolAbs::AbstractFloat, xtolAbs::AbstractFloat,
        d_inv::Vector, RSS::Float64, ntot::Real, ncoef::Int64, nsp::Int64,
        model_within::Union{Nothing, WithinSpeciesCTM}=nothing)

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
## Ancestral State Reconstruction
###############################################################################
"""
    ReconstructedStates

Type containing the inferred information about the law of the ancestral states
given the observed tips values. The missing tips are considered as ancestral states.

Reconstructed states and prediction intervals can be recovered with function `predict`,
and the standard error can be obtained with `stderror`.

The `ReconstructedStates` object has fields: `traits_nodes`, `variances_nodes`, `nodenumbers`, `traits_tips`, `tipnumbers`, `model`.
Type in "?ReconstructedStates.field" to get help on a specific field.
"""
struct ReconstructedStates
    "traits_nodes: the infered expectation of 'missing' values (ancestral nodes and missing tips)"
    traits_nodes::Vector # Nodes are actually "missing" data (including tips)
    "variances_nodes: the variance covariance matrix between all the 'missing' nodes"
    variances_nodes::Matrix
    "nodenumbers: vector of the nodes numbers, in the same order as `traits_nodes`"
    nodenumbers::Vector{Int}
    "traits_tips: the observed traits values at the tips"
    traits_tips::Vector # Observed values at tips
    "tipnumbers: vector of tips numbers, in the same order as `traits_tips`"
    tipnumbers::Vector # Observed tips only
    "model: if not missing, the `PhyloNetworkLinearModel` used for the computations."
    model::Union{PhyloNetworkLinearModel, Missing} # if empirical, corresponding fitted object
end

"""
    predict(obj::ReconstructedStates; interval::Union{Symbol,Nothing}=nothing, level::Real=0.95, text::Bool=false, digits::Int=2, missingmark::AbstractString="*", combine::Bool=false)

Estimated reconstructed states and, if `interval=:prediction`, prediction intervals with level `level` at all internal nodes and missing tips.
If `text=true`, the prediction and intervals are formated as string for easy plotting with the `nodelabel` argument to
`plot` from package [`PhyloPlots`](https://github.com/juliaphylo/PhyloPlots.jl).
In that case,
`digits` controls the number of digits shown,
`missingmark` adds a distinctive mark to prediciton of missing tips (set to `missingmark=""` for no mark),
and if `combine=true`, the prediction and bound of the intervals are combined into a single string.
"""
function StatsAPI.predict(
    obj::ReconstructedStates;
    interval::Union{Symbol,Nothing}=nothing,
    level::Real=0.95,
    text::Bool=false,
    digits::Int=2,
    missingmark::AbstractString="*",
    combine::Bool=false
    )
    res = DataFrame(
        nodenumber = [obj.nodenumbers; obj.tipnumbers], 
        prediction = [obj.traits_nodes; obj.traits_tips]
        )
    if interval === nothing
        if text 
            res[!,:prediction] = string.(round.(res[!,:prediction], digits=digits), getmissingtipmarks(obj, missingmark))
            return res
        end
        return res
    end
    if interval === :prediction
        pp = getpredint(obj; level)
        res[!,:lower] = pp[:,1]
        res[!,:upper] = pp[:,2]
        if text 
            res[!,:interval] = formatinterval(obj, res, combine, digits)
            res[!,:prediction] = string.(round.(res[!,:prediction], digits=digits), getmissingtipmarks(obj, missingmark))
            return res
        end    
        return res
    end
    throw(ArgumentError("`interval` must be one of `nothing` or `:prediction`."))
end

"""
    getmissingtipmarks(obj::ReconstructedStates, missingmark::AbstractString="*")

Create a vector of string, with a `missingmark` for tips that are missing.
"""
function getmissingtipmarks(obj::ReconstructedStates, missingmark::AbstractString="*")
    nodenumber = [obj.nodenumbers; obj.tipnumbers]
    ismissingtip = fill("", length(nodenumber))
    if !ismissing(obj.model)
        nonmissing = obj.model.nonmissing
        ind = obj.model.ind
        tipnumbers = obj.model.V.tipnumbers # all tips, even those absent from dataframe
        tipnumbers_data = tipnumbers[ind][nonmissing] # listed and data non-missing
        tipnumbers_imputed = setdiff(tipnumbers, tipnumbers_data)
        indexMissing = indexin(tipnumbers_imputed, nodenumber)
        ismissingtip[indexMissing] .*= missingmark
    end
    return(ismissingtip)
end

StatsBase.stderror(obj::ReconstructedStates) = sqrt.(diag(obj.variances_nodes))

"""
    getpredint(obj::ReconstructedStates; level::Real=0.95)

Prediction intervals with level `level` for internal nodes and missing tips.
"""
function getpredint(obj::ReconstructedStates; level::Real=0.95)
    if ismissing(obj.model)
        qq = quantile(Normal(), (1. - level)/2.)
    else
        qq = quantile(GLM.TDist(dof_residual(obj.model)), (1. - level)/2.) # TDist from Distributions
        # @warn "As the variance is estimated, the predictions intervals are not exact, and should probably be larger."
    end
    tmpnode = hcat(obj.traits_nodes, obj.traits_nodes) .+ (stderror(obj) * qq) .* [1. -1.]
    return vcat(tmpnode, hcat(obj.traits_tips, obj.traits_tips))
end

function Base.show(io::IO, obj::ReconstructedStates)
    println(io, "$(typeof(obj)):\n",
            CoefTable(hcat(vcat(obj.nodenumbers, obj.tipnumbers), vcat(obj.traits_nodes, obj.traits_tips), getpredint(obj)),
                      ["Node index", "Pred.", "Min.", "Max. (95%)"],
                      fill("", length(obj.nodenumbers)+length(obj.tipnumbers))))
end

"""
    formatinterval(obj::ReconstructedStates, pred::DataFrame, withexpectation::Bool=false, digits::Int=2)

Format the prediction intervals for the plotting function.
If `withexpectation` is set to true, then the best
predicted value is also shown along with the interval.
"""
function formatinterval(obj::ReconstructedStates, pred::DataFrame, withexpectation::Bool=false, digits::Int=2)
    pritxt = Array{AbstractString}(undef, size(pred, 1))
    for i in 1:length(obj.nodenumbers)
        !withexpectation ? sep = ", " : sep = "; " * string(round(pred[i,:prediction], digits=digits) )* "; "
        pritxt[i] = "[" * string(round(pred[i,:lower], digits=digits)) * sep * string(round(pred[i,:upper], digits=digits)) * "]"
    end
    for i in (length(obj.nodenumbers)+1):size(pred, 1)
        pritxt[i] = string(round(pred[i,:prediction], digits=digits))
    end
    return pritxt
end

#= ----- roadmap of ancestralreconstruction, continuous traits ------

all methods return a ReconstructedStates object.
core method called by every other method:

1. ancestralreconstruction(Vzz, VzyVyinvchol, RL, Y, m_y, m_z,
                                nodenumbers, tipnumbers, sigma2, add_var, model)

higher-level methods, for real data:

2. ancestralreconstruction(dataframe, net; tipnames=:tipnames, kwargs...)
   - dataframe: 2 columns only, species names & tip response values
   - fits an intercept-only model, then calls #3
   - by default without kwargs: model = BM w/o within-species variation

3. ancestralreconstruction(PhyloNetworkLinearModel[, Matrix])
   - takes a model already fitted
   - if no matrix given: the model must be intercept-only. An expanded intercept
     column is created with length = # nodes with *no* data
   - matrix: if given, must have same # of columns as the model matrix, and
     must contain the predictor(s) at nodes with *no* data, with nodes listed in
     the following order:
     * internal nodes first, in the same order in which they appear in net.node,
       i.e in V.internalnodenumbers
     * then leaves with no data, in the same order in which they appear in
       tiplabels(net), i.e. in V.tipnumbers.
   - extracts the predicted values for all network nodes, and the unscaled
     3 covariance matrices of interest (nodes with data, nodes w/o, crossed)
   - computes "universal" kriging (as opposed to "simple" kriging, which would
     simply plug-in estimates into the prediction variance formula): a term is
     added to the prediction variance, to account for the estimation of β.

methods based on simulations with a ParamsProcess "params":

4. ancestralreconstruction(net, Y, params) which calls:
   ancestralreconstruction(V::MatrixTopologicalOrder, Y, params)
   - intercept-only: known β and known variance: "simple" kriging is correct
   - BM only: params must be of type ParamsBM
=#

"""
    ancestralreconstruction(net::HybridNetwork, Y::Vector, params::ParamsBM)

Compute the conditional expectations and variances of the ancestral (un-observed)
traits values at the internal nodes of the phylogenetic network (`net`),
given the values of the traits at the tips of the network (`Y`) and some
known parameters of the process used for trait evolution (`params`, only BM with fixed root
works for now).

This function assumes that the parameters of the process are known. For a more general
function, see `ancestralreconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])`.

"""
function ancestralreconstruction(
    net::HybridNetwork,
    Y::Vector,
    params::ParamsBM
)
    V = sharedpathmatrix(net)
    ancestralreconstruction(V, Y, params)
end

function ancestralreconstruction(
    V::MatrixTopologicalOrder,
    Y::Vector,
    params::ParamsBM
)
    # Variances matrices
    Vy = V[:tips]
    Vz = V[:internalnodes]
    Vyz = V[:tipsnodes]
    R = cholesky(Vy)
    RL = R.L
    VzyVyinvchol = (RL \ Vyz)'
    # Vectors of means
    m_y = ones(size(Vy)[1]) .* params.mu # !! correct only if no predictor.
    m_z = ones(size(Vz)[1]) .* params.mu # !! works if BM no shift.
    return ancestralreconstruction(Vz, VzyVyinvchol, RL,
        Y, m_y, m_z,
        V.internalnodenumbers,
        V.tipnumbers,
        params.sigma2)
end

# Reconstruction from all the needed quantities
function ancestralreconstruction(
    Vz::Matrix,
    VzyVyinvchol::AbstractMatrix,
    RL::LowerTriangular,
    Y::Vector,
    m_y::Vector,
    m_z::Vector,
    nodenumbers::Vector,
    tipnumbers::Vector,
    sigma2::Real,
    add_var::Matrix=zeros(size(Vz)), # Additional variance for BLUP
    model::Union{PhyloNetworkLinearModel,Missing}=missing
)
    # E[z∣y] = E[z∣X] + Cov(z,y)⋅Var(y)⁻¹⋅(y-E[y∣X])
    m_z_cond_y = m_z + VzyVyinvchol * (RL \ (Y - m_y))
    V_z_cond_y = sigma2 .* (Vz - VzyVyinvchol * VzyVyinvchol')
    if !ismissing(model) && !isnothing(model.model_within) # y = last part of z
        Y = similar(Y, 0) # empty vector of similar type as Y
    end
    ReconstructedStates(m_z_cond_y, V_z_cond_y + add_var, nodenumbers, Y, tipnumbers, model)
end

#= from a fitted object: see high-level docstring below
X_n: matrix with as many columns as the number of predictors used,
     and as many rows as the number of unknown nodes or tips.

TO DO: Handle the order of internal nodes and no-data tips for matrix X_n
=#
function ancestralreconstruction(obj::PhyloNetworkLinearModel, X_n::Matrix)
    if size(X_n)[2] != length(coef(obj))
        error("""The number of predictors for the ancestral states (number of columns of X_n)
              does not match the number of predictors at the tips.""")
    end
    if size(X_n)[1] != length(obj.V.internalnodenumbers) + length(obj.V.tipnumbers)-length(obj.ind) + sum(.!obj.nonmissing)
        error("""The number of lines of the predictors does not match
              the number of nodes plus the number of missing tips.""")
    end
    #= y: observed species means at some tips
       z: trait (true species mean) at nodes to be predicted:
          - at nodes without data, i.e. internal nodes & no-data tips
          - at tips with data if within-species variation: y=ytrue+ϵ
       Vy,y = Vy,ytrue = Vytrue,y and Vytrue,z = Vyz
    =#
    m_y = predict(obj)
    m_z = X_n * coef(obj)
    # If the tips were re-organized, do the same for Vyz
    if obj.ind == [0]
        @warn """There were no indication for the position of the tips on the network.
             I am assuming that they are given in the same order.
             Please check that this is what you intended."""
        ind = collect(1:length(obj.V.tipnumbers))
    else
        ind = obj.ind
    end
    # Vyz: sharedpath. rows y: tips w/ data. cols z: internal nodes & tips w/o data
    Vyz = obj.V[:tipsnodes, ind, obj.nonmissing]
    Vzz = obj.V[:internalnodes, ind, obj.nonmissing]
    nmTipNumbers = obj.V.tipnumbers[ind][obj.nonmissing] # tips w/ data
    # no-data node numbers: for nodes (internal or tips) with no data
    ndNodeNumbers = [obj.V.internalnodenumbers; setdiff(obj.V.tipnumbers, nmTipNumbers)]
    if !isnothing(obj.model_within) # add tips with data to z
        Vtips = obj.V[:tips, ind, obj.nonmissing]
        Vzz = [Vzz Vyz'; Vyz Vtips]
        Vyz = [Vyz Vtips]
        append!(m_z, m_y)
        append!(ndNodeNumbers, nmTipNumbers)
        empty!(nmTipNumbers)
        X_n = vcat(X_n, obj.X)
    end
    VzyVyinvchol = (obj.RL \ Vyz)'
    # add_var = zeros corresponds to "simple" kriging: E[Y∣X]=Xβ with known β & variance components
    # below: "universal" kriging: β estimated, variance components known
    U = X_n - VzyVyinvchol * (obj.RL \ obj.X)
    add_var = U * vcov(obj) * U'
    @warn """These prediction intervals show uncertainty in ancestral values,
         assuming that the estimated variance rate of evolution is correct.
         Additional uncertainty in the estimation of this variance rate is
         ignored, so prediction intervals should be larger."""
    return ancestralreconstruction(
        Vzz,
        VzyVyinvchol,
        obj.RL,
        obj.Y,
        m_y,
        m_z,
        ndNodeNumbers,
        nmTipNumbers,
        sigma2_phylo(obj),
        add_var,
        obj)
end

@doc raw"""
    ancestralreconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])

Function to find the ancestral traits reconstruction on a network, given an
object fitted by function [`phylolm`](@ref). By default, the function assumes
that the regressor is just an intercept. If the value of the regressor for
all the ancestral states is known, it can be entered in X_n, a matrix with as
many columns as the number of predictors used, and as many lines as the number
of unknown nodes or tips.

Returns an object of type [`ReconstructedStates`](@ref).

# Examples

```jldoctest; filter = [r" PhyloTraits .*:\d+", ]
julia> using DataFrames, CSV # to read data file

julia> phy = readnewick(joinpath(dirname(pathof(PhyloTraits)), "..", "examples", "carnivores_tree.txt"));

julia> dat = CSV.read(joinpath(dirname(pathof(PhyloTraits)), "..", "examples", "carnivores_trait.txt"), DataFrame);

julia> using StatsModels # for statistical model formulas

julia> fitBM = phylolm(@formula(trait ~ 1), dat, phy);

julia> ancStates = ancestralreconstruction(fitBM) # Should produce a warning, as variance is unknown.
┌ Warning: These prediction intervals show uncertainty in ancestral values,
│ assuming that the estimated variance rate of evolution is correct.
│ Additional uncertainty in the estimation of this variance rate is
│ ignored, so prediction intervals should be larger.
└ @ PhyloTraits ~/build/JuliaPhylo/PhyloTraits.jl/src/traits_continuous.jl:2601
ReconstructedStates:
───────────────────────────────────────────────
  Node index      Pred.        Min.  Max. (95%)
───────────────────────────────────────────────
        -5.0   1.32139   -0.33824      2.98102
        -8.0   1.03258   -0.589695     2.65485
        -7.0   1.41575   -0.140705     2.97221
        -6.0   1.39417   -0.107433     2.89577
        -4.0   1.39961   -0.102501     2.90171
        -3.0   1.51341   -0.220523     3.24733
       -13.0   5.3192     3.92279      6.71561
       -12.0   4.51176    2.89222      6.13131
       -16.0   1.50947   -0.0186118    3.03755
       -15.0   1.67425    0.196069     3.15242
       -14.0   1.80309    0.309992     3.29618
       -11.0   2.7351     1.17608      4.29412
       -10.0   2.73217    1.12361      4.34073
        -9.0   2.41132    0.603932     4.21871
        -2.0   2.04138   -0.0340955    4.11686
        14.0   1.64289    1.64289      1.64289
         8.0   1.67724    1.67724      1.67724
         5.0   0.331568   0.331568     0.331568
         2.0   2.27395    2.27395      2.27395
         4.0   0.275237   0.275237     0.275237
         6.0   3.39094    3.39094      3.39094
        13.0   0.355799   0.355799     0.355799
        15.0   0.542565   0.542565     0.542565
         7.0   0.773436   0.773436     0.773436
        10.0   6.94985    6.94985      6.94985
        11.0   4.78323    4.78323      4.78323
        12.0   5.33016    5.33016      5.33016
         1.0  -0.122604  -0.122604    -0.122604
        16.0   0.73989    0.73989      0.73989
         9.0   4.84236    4.84236      4.84236
         3.0   1.0695     1.0695       1.0695
───────────────────────────────────────────────

julia> using StatsBase # for predict function

julia> predict(ancStates)
31×2 DataFrame
 Row │ nodenumber  prediction 
     │ Int64       Float64    
─────┼────────────────────────
   1 │         -5    1.32139
   2 │         -8    1.03258
   3 │         -7    1.41575
   4 │         -6    1.39417
   5 │         -4    1.39961
   6 │         -3    1.51341
   7 │        -13    5.3192
   8 │        -12    4.51176
  ⋮  │     ⋮           ⋮
  25 │         10    6.94985
  26 │         11    4.78323
  27 │         12    5.33016
  28 │          1   -0.122604
  29 │         16    0.73989
  30 │          9    4.84236
  31 │          3    1.0695
               16 rows omitted

julia> predict(ancStates, interval = :prediction)
31×4 DataFrame
 Row │ nodenumber  prediction  lower       upper     
     │ Int64       Float64     Float64     Float64   
─────┼───────────────────────────────────────────────
   1 │         -5    1.32139   -0.33824     2.98102
   2 │         -8    1.03258   -0.589695    2.65485
   3 │         -7    1.41575   -0.140705    2.97221
   4 │         -6    1.39417   -0.107433    2.89577
   5 │         -4    1.39961   -0.102501    2.90171
   6 │         -3    1.51341   -0.220523    3.24733
   7 │        -13    5.3192     3.92279     6.71561
   8 │        -12    4.51176    2.89222     6.13131
  ⋮  │     ⋮           ⋮           ⋮           ⋮
  25 │         10    6.94985    6.94985     6.94985
  26 │         11    4.78323    4.78323     4.78323
  27 │         12    5.33016    5.33016     5.33016
  28 │          1   -0.122604  -0.122604   -0.122604
  29 │         16    0.73989    0.73989     0.73989
  30 │          9    4.84236    4.84236     4.84236
  31 │          3    1.0695     1.0695      1.0695
                                      16 rows omitted

julia> using PhyloPlots # next: plot ancestral states on the tree

julia> plot(phy, nodelabel = predict(ancStates));

julia> pred = predict(ancStates, interval = :prediction, text = true);

julia> plot(phy, nodelabel = pred[!,[:nodenumber,:interval]]);

julia> allowmissing!(dat, :trait);

julia> dat[[2, 5], :trait] .= missing; # missing values allowed to fit model

julia> fitBM = phylolm(@formula(trait ~ 1), dat, phy);

julia> ancStates = ancestralreconstruction(fitBM);
┌ Warning: These prediction intervals show uncertainty in ancestral values,
│ assuming that the estimated variance rate of evolution is correct.
│ Additional uncertainty in the estimation of this variance rate is
│ ignored, so prediction intervals should be larger.
└ @ PhyloTraits ~/build/JuliaPhylo/PhyloTraits.jl/src/traits_continuous.jl:2601

julia> first(predict(ancStates), 3) # looking at first 3 nodes only
3×2 DataFrame
 Row │ nodenumber  prediction 
     │ Int64       Float64    
─────┼────────────────────────
   1 │         -5     1.42724
   2 │         -8     1.35185
   3 │         -7     1.61993

julia> first(predict(ancStates, interval=:prediction), 3)
3×4 DataFrame
 Row │ nodenumber  prediction  lower      upper   
     │ Int64       Float64     Float64    Float64 
─────┼────────────────────────────────────────────
   1 │         -5     1.42724  -0.373749  3.22824
   2 │         -8     1.35185  -0.698432  3.40214
   3 │         -7     1.61993  -0.17179   3.41165

julia> plot(phy, nodelabel = predict(ancStates, text=true));

julia> pred = predict(ancStates, interval = :prediction, text = true);

julia> plot(phy, nodelabel = pred[!,[:nodenumber,:interval]]);
```
"""
function ancestralreconstruction(obj::PhyloNetworkLinearModel)
    # default reconstruction for known predictors
    if ((size(obj.X)[2] != 1) || !any(obj.X .== 1)) # Test if the regressor is just an intercept.
        error("""Predictor(s) other than a plain intercept are used in this `PhyloNetworkLinearModel` object.
    These predictors are unobserved at ancestral nodes, so they cannot be used
    for the ancestral state reconstruction. If these ancestral predictor values
    are known, please provide them as a matrix argument to the function.
    Otherwise, you might consider doing a multivariate linear regression (not implemented yet).""")
    end
    X_n = ones((length(obj.V.nodenumbers_toporder) - sum(obj.nonmissing), 1))
    ancestralreconstruction(obj, X_n)
end

"""
    ancestralreconstruction(fr::AbstractDataFrame, net::HybridNetwork; kwargs...)

Estimate the ancestral traits on a network, given some data at the tips.
Uses function [`phylolm`](@ref) to perform a phylogenetic regression of the data against an
intercept (amounts to fitting an evolutionary model on the network).

See documentation on [`phylolm`](@ref) and `ancestralreconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix])`
for further details.

Returns an object of type [`ReconstructedStates`](@ref).
"""
function ancestralreconstruction(
    fr::AbstractDataFrame,
    net::HybridNetwork;
    tipnames::Symbol=:tipnames,
    kwargs...
)
    nn = DataFrames.propertynames(fr)
    datpos = nn .!= tipnames
    if sum(datpos) > 1
        error("""Besides one column labelled '$tipnames', the dataframe fr should have
              only one column, corresponding to the data at the tips of the network.""")
    end
    f = @eval(@formula($(nn[datpos][1]) ~ 1))
    reg = phylolm(f, fr, net; tipnames=tipnames, kwargs...)
    return ancestralreconstruction(reg)
end






#################################################
## Old version of phylolm (naive)
#################################################

# function phylolmNaive(X::Matrix, Y::Vector, net::HybridNetwork, model::AbstractString="BM")
#   # Geting variance covariance
#   V = sharedpathmatrix(net)
#   Vy = extractVarianceTips(V, net)
#   # Needed quantities (naive)
#   ntaxa = length(Y)
#   Vyinv = inv(Vy)
#   XtVyinv = X' * Vyinv
#   logdetVy = logdet(Vy)
#        # beta hat
#   betahat = inv(XtVyinv * X) * XtVyinv * Y
#        # sigma2 hat
#   fittedValues =  X * betahat
#   residuals = Y - fittedValues
#   sigma2hat = 1/ntaxa * (residuals' * Vyinv * residuals)
#        # log likelihood
#   loglik = - 1 / 2 * (ntaxa + ntaxa * log(2 * pi) + ntaxa * log(sigma2hat) + logdetVy)
#   # Result
# # res = phyloNetworkRegression(betahat, sigma2hat[1], loglik[1], V, Vy, fittedValues, residuals)
#   return((betahat, sigma2hat[1], loglik[1], V, Vy, logdetVy, fittedValues, residuals))
# end
