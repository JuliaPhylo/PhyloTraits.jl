"""
    posterior_logtreeweight(obj::SSM, trait = 1)

Array A of log-posterior probabilities for each tree displayed in the network:
such that A[t] = log of P(tree `t` | trait `trait`)
if a single `trait` is requested, or A[t,i]= log of P(tree `t` | trait `i`)
if `trait` is a vector or range (e.g. `trait = 1:obj.nsites`).
These probabilities are conditional on the model parameters in `obj`.

Displayed trees are listed in the order in which they are stored in the fitted
model object `obj`.

**Precondition**: `_loglikcache` updated by [`discrete_corelikelihood!`](@ref)

# examples

```jldoctest
julia> net = readnewick("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]); # arbitrary rates

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit = fitdiscrete(net, m1, dat); # optimized rates: α=0.27 and β=0.35

julia> pltw = PhyloTraits.posterior_logtreeweight(fit);

julia> round.(exp.(pltw), digits=5) # posterior trees probabilities (sum up to 1)
2-element Vector{Float64}:
 0.91983
 0.08017

julia> round.(exp.(fit.priorltw), digits=4) # the prior tree probabilities are similar here (tiny data set!)
2-element Vector{Float64}:
 0.9
 0.1
```
"""
function posterior_logtreeweight(obj::SSM, trait = 1)
    # ts[site,tree] = log P(data and tree) at site, integrated over rates
    d = length(size(trait)) # 0 if single trait, 1 if vector of several traits
    ts = dropdims(mapslices(logsumexp, view(obj._loglikcache, trait,:,:),
                            dims=d+1); dims=1)
    if d>0 ts = permutedims(ts); end # now: ts[tree] or ts[tree,site]
    siteliks = mapslices(logsumexp, ts, dims=1) # 1 x ntraits array (or 1-element vector)
    ts .-= siteliks
    return ts
end

"""
    posterior_loghybridweight(obj::SSM, hybrid_name, trait = 1)
    posterior_loghybridweight(obj::SSM, edge_number, trait = 1)

Log-posterior probability for all trees displaying the minor parent edge
of hybrid node named `hybrid_name`, or displaying the edge number `edge_number`.
That is: log of P(hybrid minor parent | trait) if a single `trait` is requested,
or A[i]= log of P(hybrid minor parent | trait `i`)
if `trait` is a vector or range (e.g. `trait = 1:obj.nsites`).
These probabilities are conditional on the model parameters in `obj`.

**Precondition**: `_loglikcache` updated by [`discrete_corelikelihood!`](@ref)

# examples

```jldoctest
julia> net = readnewick("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]); # arbitrary rates

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit = fitdiscrete(net, m1, dat); # optimized rates: α=0.27 and β=0.35

julia> plhw = PhyloTraits.posterior_loghybridweight(fit, "H1");

julia> round(exp(plhw), digits=5) # posterior probability of going through minor hybrid edge
0.08017

julia> hn = net.node[3]; getparentedgeminor(hn).gamma # prior probability
0.1
```
"""
function posterior_loghybridweight(obj::SSM, hybridname::String, trait = 1)
    hn_index = findfirst(n -> n.name == hybridname, obj.net.node)
    isnothing(hn_index) && error("node named $hybridname not found")
    hn = obj.net.node[hn_index]
    hn.hybrid || error("node named $hybridname is not a hybrid node")
    me = getparentedgeminor(hn)
    posterior_loghybridweight(obj, me.number, trait)
end
function posterior_loghybridweight(obj::SSM, edgenum::Integer, trait = 1)
    tpp = posterior_logtreeweight(obj, trait) # size: (ntree,) or (ntree,ntraits)
    hasedge = tree -> any(e.number == edgenum for e in tree.edge)
    tokeep = map(hasedge, obj.displayedtree)
    tppe = view(tpp, tokeep, :) # makes it a matrix
    epp = dropdims(mapslices(logsumexp, tppe, dims=1); dims=2)
    return (size(epp)==(1,) ? epp[1] : epp) # scalar or vector
end

"""
    ancestralreconstruction(obj::SSM, trait::Integer = 1)

Estimate the marginal probability of ancestral states for discrete character
number `trait` (first trait by default).
The parameters of the [`StatisticalSubstitutionModel`](@ref) object `obj`
must first be fitted using [`fitdiscrete`](@ref), and ancestral state reconstruction
is conditional on the estimated parameters. If these parameters were estimated
using all traits, they are used as is, to do ancestral state reconstruction of the
particular `trait` of interest.

**output**: data frame with a first column for the node numbers, a second column for
the node labels, and a column for each possible state: the entries in these columns
give the marginal probability that a given node has a given state.

warnings:
- node numbers and node labels refer to those in `obj.net`, which might
  have a different internal representation of nodes than the original network
  used to build `obj`.
- `obj` is modified: its likelihood fields (forward, directional & backward)
  are updated to make sure that they correspond to the current parameter values
  in `obj.model`, and to the `trait` of interest.

limitations: the following are not checked.
- Assumes that every node in the large network is also present
  (with descendant leaves) in each displayed tree.
  This is not true if the network is not tree-child...
- Assumes that the root is also in each displayed tree, which
  may not be the case if the root had a hybrid child edge.

See also [`posterior_logtreeweight`](@ref) and
[`discrete_backwardlikelihood_trait!`](@ref) to update `obj.backwardlik`.

# examples

```jldoctest
julia> net = readnewick("(((A:2.0,(B:1.0)#H1:0.1::0.9):1.5,(C:0.6,#H1:1.0::0.1):1.0):0.5,D:2.0);");

julia> m1 = BinaryTraitSubstitutionModel([0.1, 0.1], ["lo", "hi"]);

julia> using DataFrames

julia> dat = DataFrame(species=["C","A","B","D"], trait=["hi","lo","lo","hi"]);

julia> fit1 = fitdiscrete(net, m1, dat);

julia> asr = ancestralreconstruction(fit1)
9×4 DataFrame
 Row │ nodenumber  nodelabel  lo        hi
     │ Int64       String     Float64   Float64
─────┼───────────────────────────────────────────
   1 │          1  A          1.0       0.0
   2 │          2  B          1.0       0.0
   3 │          3  C          0.0       1.0
   4 │          4  D          0.0       1.0
   5 │          5  5          0.286021  0.713979
   6 │          6  6          0.319456  0.680544
   7 │          7  7          0.16855   0.83145
   8 │          8  8          0.767359  0.232641
   9 │          9  H1         0.782776  0.217224

julia> using PhyloPlots

julia> plot(fit1.net, nodelabel = asr[!,[:nodenumber, :lo]], tipoffset=0.2); # pp for "lo" state
```
"""
function ancestralreconstruction(obj::SSM, trait::Integer = 1)
    # posterior probability of state i at node n: proportional to
    # sum_{tree t, rate r} exp( ltw[t] + backwardll[i,n] given t,r + forwardll[i,n] given t,r ) / nr
    trait <= obj.nsites || error("trait $trait is larger than the number of traits in the data")
    nnodes = length(obj.net.node) # may be smaller than 2nd size of bkd or frd
    update_logtrans(obj)
    bkd = view(obj.backwardlik, :, 1:nnodes)
    frd = view(obj.forwardlik, :, 1:nnodes)
    ltw = obj.priorltw
    res = similar(bkd) # first: hold the cumulative logsumexp of bkd + frd + ltw
    fill!(res, -Inf64)
    nr = length(obj.ratemodel.ratemultiplier)
    lrw = obj.ratemodel.lograteweight
    for t in 1:length(obj.displayedtree)
        ltprior = ltw[t]
        for ri in 1:nr
            # update forward & directional likelihoods
            discrete_corelikelihood_trait!(obj,t,trait,ri)
            # update backward likelihoods
            discrete_backwardlikelihood_trait!(obj,t,ri)
            # P{state i at node n} ∝ bkd[i,n] * frd[i,n] given tree & rate:
            # res = logaddexp(res, ltw[t] + lrw[ri] + bkd + frd)
            broadcast!(logaddexp, res, res, (ltprior + lrw[ri]) .+ bkd .+ frd)
        end
    end
    # normalize the results at each node: p_i / sum(p_j over all states j)
    traitloglik = logsumexp(res[:,1]) # sum p_j at node 1 or at any node = loglikelihood
    res .= exp.(res .- traitloglik)
    nodestringlabels = Vector{String}(undef, nnodes)
    for n in obj.net.node
        nodestringlabels[n.number] = (n.name == "" ? string(n.number) : n.name)
    end
    dat = DataFrame(transpose(res), Symbol.(getlabels(obj.model)))
    insertcols!(dat, 1, :nodenumber => collect(1:nnodes), makeunique=true)
    insertcols!(dat, 2, :nodelabel  => nodestringlabels,  makeunique=true)
    return dat
end

"""
    discrete_backwardlikelihood_trait!(obj::SSM, tree::Integer, ri::Integer)

Update and return the backward likelihood (last argument `backwardlik`)
assuming rate category `ri` and tree index `tree`,
using current forward and backwards likelihoods in `obj`:
these depend on the trait (or site) given to the last call to
`discrete_corelikelihood_trait!`.
Used by `ancestralreconstruction`.

**warning**: assume correct transition probabilities.
"""
function discrete_backwardlikelihood_trait!(obj::SSM, t::Integer, ri::Integer)
    backwardlik = obj.backwardlik
    directlik  = obj.directlik
    tree = obj.displayedtree[t]
    k = nstates(obj.model)
    fill!(backwardlik, 0.0) # re-initialize for each trait, each iteration
    bkwtmp = Vector{Float64}(undef, k) # to hold bkw lik without parent edge transition
    if typeof(obj.model) <: NASM
        logprior = log.(stationary(obj.model))
    else #trait models
        logprior = [-log(k) for i in 1:k] # uniform prior at root
    end
    for ni in 1:length(tree.vec_node) # pre-order traversal to calculate backwardlik
        n = tree.vec_node[ni]
        nnum = n.number
        if ni == 1 # n is the root
            backwardlik[:,nnum] = logprior
        else
            pe = getparentedge(n)
            pn = getparent(pe)
            bkwtmp[:] = backwardlik[:,pn.number] # use bktmp's original memory
            for se in pn.edge
                if se != pe && pn == getparent(se) # then se is sister edge to pe
                    bkwtmp .+= view(directlik, :,se.number)
                end
            end
            lt = view(obj.logtrans, :,:,pe.number,ri)
            for j in 1:k # state at node n
                backwardlik[j,nnum] = logsumexp(bkwtmp + view(lt,:,j))
            end
        end
    end
    return backwardlik
end
