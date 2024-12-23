"""
    rand([rng::AbstractRNG,]
          model::TraitSubstitutionModel,
          t::Float64,
          start::AbstractVector{Int})

Simulate discrete traits along one edge of length t. A random number generator
`rng` is optional. `start` must be a vector of integers, each representing the
starting value of one trait.

# examples

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0)
Binary Trait Substitution Model:
rate 0→1 α=1.0
rate 1→0 β=2.0

julia> using StableRNGs; rng = StableRNG(135);

julia> rand(rng, m1, 0.2, [1,2,1,2,1])
5-element Vector{Int64}:
 2
 2
 1
 2
 1
```
"""
function rand(obj::SubstitutionModel, t::Float64, start::AbstractVector{Int})
    rand(default_rng(), obj, t, start)
end
function rand(
    rng::AbstractRNG,
    obj::SubstitutionModel,
    t::Float64,
    start::AbstractVector{Int}
)
    res = Vector{Int}(undef, length(start))
    rand!(rng, res, obj, t, start)
end

"""
    rand!(rng::AbstractRNG,
          end::AbstractVector{Int},
          model::TraitSubstitutionModel,
          t::Float64,
          start::AbstractVector{Int})

Simulate discrete traits along one edge of length `t`
like [`rand`](@ref PhyloTraits.rand),
but modifying `end` in place to store the simulated values.
"""
function rand!(
    rng::AbstractRNG,
    endTrait::AbstractVector{Int},
    obj::SubstitutionModel,
    t::Float64,
    start::AbstractVector{Int}
)
    Pt = P(obj, t)
    k = size(Pt, 1) # number of states
    w = [aweights(Pt[i,:]) for i in 1:k]
    for i in eachindex(start)
        endTrait[i] = sample(rng, 1:k, w[start[i]])
    end
    return endTrait
end

"""
    rand([rng::AbstractRNG,]
          model::TraitSubstitutionModel,
          net::HybridNetwork;
          ntraits=1,
          keepinternal=true,
          checkpreorder=true)

Simulate evolution of discrete traits on a rooted evolutionary network based on
the supplied evolutionary model. Trait sampling is uniform at the root.

optional arguments:

- `ntraits`: number of traits to be simulated (default: 1 trait).
- `keepinternal`: if true, export character states at all nodes, including
  internal nodes. if false, export character states at tips only.

output:

- matrix of character states with one row per trait, one column per node;
  these states are *indices* in `model.label`, not the trait labels themselves.
- vector of node labels (for tips) or node numbers (for internal nodes)
  in the same order as columns in the character state matrix

# examples

```jldoctest
julia> m1 = BinaryTraitSubstitutionModel(1.0, 2.0, ["low","high"]);

julia> net = readnewick("(((A:4.0,(B:1.0)#H1:1.1::0.9):0.5,(C:0.6,#H1:1.0::0.1):1.0):3.0,D:5.0);");

julia> using Random; Random.seed!(95);

julia> trait, lab = rand(m1, net)
([1 2 … 1 1], ["-2", "D", "-3", "-6", "C", "-4", "H1", "B", "A"])

julia> trait
1×9 Matrix{Int64}:
 1  2  1  1  2  2  1  1  1

julia> lab
9-element Vector{String}:
 "-2" 
 "D"  
 "-3" 
 "-6" 
 "C"  
 "-4" 
 "H1"
 "B"  
 "A"  
```
"""
function rand(obj::SubstitutionModel, net::HybridNetwork; kwargs...)
    rand(default_rng(), obj, net; kwargs...)
end
function rand(
    rng::AbstractRNG,
    obj::SubstitutionModel,
    net::HybridNetwork;
    ntraits::Int=1,
    keepinternal::Bool=true,
    checkpreorder::Bool=true
)
    net.isrooted || error("net needs to be rooted for preorder recursion")
    checkpreorder && preorder!(net)
    nnodes = net.numnodes
    M = Matrix{Int}(undef, ntraits, nnodes) # M[i,j]= trait i for node j
    rand!(rng, M, obj, net)
    if !keepinternal
        M = PN.getTipSubmatrix(M, net, indexation=:cols) # subset columns only. rows=traits
        nodeLabels = [n.name for n in net.vec_node if n.leaf]
    else
        nodeLabels = [n.name == "" ? string(n.number) : n.name for n in net.vec_node]
    end
    return M, nodeLabels
end

function rand!(
    rng::AbstractRNG,
    M::Matrix,
    obj::SubstitutionModel,
    net::HybridNetwork
)
    return PN.traversal_preorder!(
        net.vec_node,
        M, # updates M in place
        updateRootRandomTrait!,
        updateTreeRandomTrait!,
        updateHybridRandomTrait!,
        obj,
        rng)
end

function updateRootRandomTrait!(V::AbstractArray, i::Int, obj, rng)
    sample!(rng, 1:nstates(obj), view(V, :, i)) # uniform at the root
    return true
end

function updateTreeRandomTrait!(
    V::Matrix,
    i::Int,
    parentIndex::Int,
    edge::Edge,
    obj,
    rng,
)
    rand!(rng, view(V, :, i), obj, edge.length, view(V, :, parentIndex))
    return true
end

function updateHybridRandomTrait!(
    V::Matrix,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    obj,
    rng,
)
    nump = length(parindx) # 2 parents if bicombining
    rand!(rng, view(V, :, i), obj, paredge[1].length, view(V, :, parindx[1]))
    tmp = [rand(rng, obj, paredge[p].length, view(V, :, parindx[p])) for p in 2:nump]
    cs = cumsum(e.gamma for e in paredge) # last value should be 1 = sum of γs
    for j in 1:size(V,1) # loop over traits
        u = rand(rng) # next: find index p such that cs[p-1] < u < cs[p]
        p = findfirst(s -> u < s, cs) # inherit from parent p
        if p > 1 # parent 1 was stored in V already: nothing to do if p=1
            V[j,i] = tmp[p-1][j] # switch to inherit trait of parent p
        end
    end
    return true
end
