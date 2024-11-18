"""
    TraitSimulation

Result of a trait simulation on an [`HybridNetwork`](@ref) with function [`simulate`](@ref).

The following functions and extractors can be applied to it: [`tiplabels`](@ref), `obj[:Tips]`, `obj[:InternalNodes]` (see documentation for function [`getindex(::TraitSimulation, ::Symbol)`](@ref)).

The `TraitSimulation` object has fields: `M`, `params`, `model`.
"""
struct TraitSimulation
    M::MatrixTopologicalOrder
    params::ParamsProcess
    evomodel::AbstractString
end

function Base.show(io::IO, obj::TraitSimulation)
    disp = "$(typeof(obj)):\n"
    disp = disp * "Trait simulation results on a network with $(length(obj.M.tipNames)) tips, using a $(obj.evomodel) model, with parameters:\n"
    disp = disp * paramstable(obj.params)
    println(io, disp)
end

tiplabels(obj::TraitSimulation) = tiplabels(obj.M)

"""
    simulate([rng::AbstractRNG,]
             net::HybridNetwork,
             params::ParamsProcess,
             checkpreorder::Bool=true)

Simulate traits on `net` using the parameters `params`. For now, only
parameters of type [`ParamsBM`](@ref) (univariate Brownian Motion) and
[`ParamsMultiBM`](@ref) (multivariate Brownian motion) are accepted.

The simulation using a recursion from the root to the tips of the network,
therefore, a pre-ordering of nodes is needed. If `checkpreorder=true` (default),
[`preorder!`](@ref) is called on the network beforehand. Otherwise, it is assumed
that the preordering has already been calculated.

Returns an object of type [`TraitSimulation`](@ref),
which has a matrix with the trait expecations and simulated trait values at
all the nodes.

See examples below for accessing expectations and simulated trait values.

# Examples

## Univariate

```jldoctest
julia> phy = readnewick("(A:2.5,((U:1,#H1:0.5::0.4):1,(C:1,(D:0.5)#H1:0.5::0.6):1):0.5);");

julia> par = ParamsBM(1, 0.1) # BM with expectation 1 and variance 0.1.
ParamsBM:
Parameters of a BM with fixed root:
mu: 1
Sigma2: 0.1


julia> using Random; Random.seed!(17920921); # for reproducibility

julia> sim = simulate(phy, par) # Simulate on the tree.
TraitSimulation:
Trait simulation results on a network with 4 tips, using a BM model, with parameters:
mu: 1
Sigma2: 0.1


julia> traits = sim[:Tips] # Extract simulated values at the tips.
4-element Vector{Float64}:
 0.9664650558470932
 0.4104321932336118
 0.2796524923704289
 0.7306692819731366

julia> sim.M.tipNames # name of tips, in the same order as values above
4-element Vector{String}:
 "A"
 "U"
 "C"
 "D"

julia> traits = sim[:InternalNodes] # Extract simulated values at internal nodes. Order: as in sim.M.internalNodeNumbers
5-element Vector{Float64}:
 0.5200361297500204
 0.8088890626285765
 0.9187604100796469
 0.711921371091375
 1.0

julia> traits = sim[:All] # simulated values at all nodes, ordered as in sim.M.nodeNumbersTopOrder
9-element Vector{Float64}:
 1.0
 0.711921371091375
 0.9187604100796469
 0.2796524923704289
 0.5200361297500204
 0.8088890626285765
 0.7306692819731366
 0.4104321932336118
 0.9664650558470932

julia> traits = sim[:Tips, :Exp] # Extract expected values at the tips (also works for sim[:All, :Exp] and sim[:InternalNodes, :Exp]).
4-element Vector{Float64}:
 1.0
 1.0
 1.0
 1.0
```

## Multivariate

```jldoctest
julia> phy = readnewick("(A:2.5,((B:1,#H1:0.5::0.4):1,(C:1,(V:0.5)#H1:0.5::0.6):1):0.5);");

julia> par = ParamsMultiBM([1.0, 2.0], [1.0 0.5; 0.5 1.0]) # BM with expectation [1.0, 2.0] and variance [1.0 0.5; 0.5 1.0].
ParamsMultiBM:
Parameters of a MBD with fixed root:
mu: [1.0, 2.0]
Sigma: [1.0 0.5; 0.5 1.0]

julia> using Random; Random.seed!(17920921); # for reproducibility

julia> sim = simulate(phy, par) # simulate on the phylogeny
TraitSimulation:
Trait simulation results on a network with 4 tips, using a MBD model, with parameters:
mu: [1.0, 2.0]
Sigma: [1.0 0.5; 0.5 1.0]


julia> traits = sim[:Tips] # Extract simulated values at the tips (each column contains the simulated traits for one node).
2×4 Matrix{Float64}:
 2.99232  -0.548734  -1.79191  -0.773613
 4.09575   0.712958   0.71848   2.00343

julia> traits = sim[:InternalNodes] # simulated values at internal nodes. order: same as in sim.M.internalNodeNumbers
2×5 Matrix{Float64}:
 -0.260794  -1.61135  -1.93202   0.0890154  1.0
  1.46998    1.28614   0.409032  1.94505    2.0

julia> traits = sim[:All]; # 2×9 Matrix: values at all nodes, ordered as in sim.M.nodeNumbersTopOrder

julia> sim[:Tips, :Exp] # Extract expected values (also works for sim[:All, :Exp] and sim[:InternalNodes, :Exp])
2×4 Matrix{Float64}:
 1.0  1.0  1.0  1.0
 2.0  2.0  2.0  2.0
```
"""
function simulate(net::HybridNetwork, params::ParamsProcess, checkpreorder::Bool=true)
    simulate(default_rng(), net, params, checkpreorder)
end
function simulate(
    rng::AbstractRNG,
    net::HybridNetwork,
    params::ParamsProcess,
    checkpreorder::Bool=true
)
    if isa(params, ParamsBM)
        model = "BM"
    elseif isa(params, ParamsMultiBM)
        model = "MBD"
    else
        error("The 'simulate' function only works for a BM process (for now).")
    end
    !ismissing(params.shift) || (params.shift = ShiftNet(net, process_dim(params)))

    net.isrooted || error("The net needs to be rooted for trait simulation.")
    !anyShiftOnRootEdge(params.shift) || error("Shifts are not allowed above the root node. Please put all root specifications in the process parameter.")

    checkpreorder && preorder!(net)
    f = preorderFunctions(params, rng)
    V = PN.traversal_preorder(net.vec_node,
            f[:init], f[:root], f[:tree], f[:hybrid], params)
    M = MatrixTopologicalOrder(V, net, :c) # nodes in columns of V
    TraitSimulation(M, params, model)
end


function preorderFunctions(::ParamsBM, rng::AbstractRNG)
    return (init = initSimulateBM,
            root = updateRootSimulateBM!(rng),
            tree = updateTreeSimulateBM!(rng),
            hybrid = updateHybridSimulateBM!(rng))
end
function preorderFunctions(::ParamsMultiBM, rng::AbstractRNG)
    return (init = initSimulateMBD,
            root = updateRootSimulateMBD!(rng),
            tree = updateTreeSimulateMBD!(rng),
            hybrid = updateHybridSimulateMBD!(rng))
end


function anyShiftOnRootEdge(shift::ShiftNet)
    nodInd = getShiftRowInds(shift)
    for n in shift.net.vec_node[nodInd]
        !(getMajorParentEdgeNumber(n) == -1) || return(true)
    end
    return false
end

# Initialization of the structure
function initSimulateBM(nodes::Vector{Node}, ::ParamsBM)
    return(zeros(2, length(nodes)))
end

function initSimulateMBD(nodes::Vector{Node}, params::ParamsMultiBM)
    p = process_dim(params)
    return zeros(2 * p, length(nodes)) # [means vals]
end

# Initialization of the root
function updateRootSimulateBM!(rng::AbstractRNG)
    f = function(M::Matrix, i::Int, params::ParamsBM)
        M[1, i] = params.mu # expectation
        M[2, i] = params.mu
        if params.randomRoot
            M[2, i] += sqrt(params.varRoot) * randn(rng) # random value
        end
        return true
    end
    return f
end
function updateRootSimulateMBD!(rng::AbstractRNG)
    f = function(M::Matrix{Float64}, i::Int, params::ParamsMultiBM)
        p = process_dim(params)
        means, vals = partitionMBDMatrix(M, p)
        means[:, i] .= params.mu # expectation
        vals[:, i] .= params.mu
        if params.randomRoot
            vals[:, i] += cholesky(params.varRoot).L * randn(rng, p) # random value
        end
        return true
    end
    return f
end

# Going down to a tree node
function updateTreeSimulateBM!(rng::AbstractRNG)
    f = function(M::Matrix,
                 i::Int,
                 parentIndex::Int,
                 edge::Edge,
                 params::ParamsBM)
        M[1, i] = M[1, parentIndex] + params.shift.shift[i] # expectation
        M[2, i] = M[2, parentIndex] + params.shift.shift[i] +
            sqrt(params.sigma2 * edge.length) * randn(rng) # random value
        return true
    end
    return f
end
function updateTreeSimulateMBD!(rng::AbstractRNG)
    f = function(M::Matrix{Float64},
                 i::Int,
                 parentIndex::Int,
                 edge::Edge,
                 params::ParamsMultiBM)
        p = process_dim(params)
        means, vals = partitionMBDMatrix(M, p)
        μ = @view means[:, i]
        val = @view vals[:, i]
        # μ .= means[:, parentIndex] + params.shift.shift[i, :]
        μ .= @view means[:, parentIndex]
        μ .+= @view params.shift.shift[i, :]
        # val .= sqrt(edge.length) * params.L * randn(rng, p) + vals[:, parentIndex] + params.shift.shift[i, :]
        mul!(val, params.L, randn(rng, p))
        val .*= sqrt(edge.length)
        val .+= @view vals[:, parentIndex]
        val .+= params.shift.shift[i, :]
        return true
    end
    return f
end

# Going down to an hybrid node
function updateHybridSimulateBM!(rng::AbstractRNG)
    f = function(M::Matrix,
                 i::Int,
                 parindx::AbstractVector{Int},
                 paredge::AbstractVector{Edge},
                 params::ParamsBM)
        iter = zip(parindx, paredge)
        M[1, i] = sum(e.gamma *  M[1,j] for (j,e) in iter) # expectation
        M[2, i] = sum(e.gamma * (M[2,j] + sqrt(params.sigma2 * e.length) * randn(rng)) for (j,e) in iter)
        return true
    end
    return f
end
# assumes a bicombining network for multivariate
function updateHybridSimulateMBD!(rng::AbstractRNG)
    f = function(M::Matrix{Float64},
                 i::Int,
                 parindx::AbstractVector{Int},
                 paredge::AbstractVector{Edge},
                 params::ParamsMultiBM)
        p = process_dim(params)
        means, vals = partitionMBDMatrix(M, p)
        μ = @view means[:, i]
        val = @view vals[:, i]
        μ1 = @view means[:, parindx[1]]
        μ2 = @view means[:, parindx[2]]
        v1 = @view vals[:, parindx[1]]
        v2 = @view vals[:, parindx[2]]
        # means[:, i] .= edge1.gamma * μ1 + edge2.gamma * μ2
        mul!(μ, μ1, paredge[1].gamma)
        BLAS.axpy!(paredge[2].gamma, μ2, μ)  # expectation
        # val .=  edge1.gamma * (v1 + sqrt(edge1.length) * params.L * r1) +
        #         edge2.gamma * (v2 + sqrt(edge2.length) * params.L * r2) # random value
        mul!(val, params.L, randn(rng, p))
        val .*= sqrt(paredge[1].length)
        val .+= v1
        buffer = params.L * randn(rng, p)
        buffer .*= sqrt(paredge[2].length)
        buffer .+= v2
        BLAS.axpby!(paredge[2].gamma, buffer, paredge[1].gamma, val) # random value
        return true
    end
    return f
end

"""
    getindex(obj, d)

Getting submatrices of an object of type [`TraitSimulation`](@ref).

# Arguments
* `obj::TraitSimulation`: the matrix from which to extract.
* `d::Symbol`: a symbol precising which sub-matrix to extract. Can be:
  * `:Tips` columns and/or rows corresponding to the tips
  * `:InternalNodes` columns and/or rows corresponding to the internal nodes
"""
function Base.getindex(obj::TraitSimulation, d::Symbol, w::Symbol=:Sim)
    inds = siminds(obj.params, w)
    return getindex(obj.M, d)[inds, :]
end

function siminds(::ParamsBM, w::Symbol)
    if w == :Sim
        return 2
    elseif w == :Exp
        return 1
    else
        error("The argument 'w' must be ':Sim' or ':Exp'. (':$w' was supplied)")
    end
end

function siminds(params::ParamsMultiBM, w::Symbol)
    p = process_dim(params)
    if w == :Sim
        return (p + 1):(2 * p)
    elseif w == :Exp
        return 1:p
    else
        error("The argument 'w' must be ':Sim' or ':Exp'. (':$w' was supplied)")
    end
end
