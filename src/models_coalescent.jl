
"""
    GaussianCoalescent(v0,σ²,Ne)

[`ContinuousTraitEM`](@ref) type to model the evolution of a polygenic trait X
with variance v0 at the root population, and variance rate σ².
The model assumes that X is controlled by a large number L of loci, each
locus l having an infinitesimally small effect Yl on the trait,
acting additively across loci: `X = (sum_l Y_l)/√L`.
The normalization by √L is to simplify notations: the effect of locus l is
`Y_l/√L`, and becomes infinitesimally small as L grows to infinity.

`v0` is the variance of each `Y_l`, and of `X` in the root population.
Each `Y_l` evolves along each own gene tree, arising from the population
phylogeny according to the network multispecies coalescent process.
Conditional on its gene tree, `Y_l` evolves according to a centered process
with increments of variance σ²t over t generations,
such as a Brownian motion process, a compound Poisson process,
or a variance-gamma process.

Ne is the haploid effective population size (that is, 2N in the standard case
of autosomal loci in diploid organisms).
Setting Ne=1 corresponds to using branch lengths in coalescent units in the
species phylogeny, and interpreting σ² as a variance rate per coalescent unit.

Ne and σ² are assumed constant across populations.
"""
struct GaussianCoalescent <: ContinuousTraitEM
    "variance at the root population"
    v0::Float64
    "variance rate σ² per generation for the polygenic trait (rate σ²/L for each of L loci)"
    sigma2::Float64 # between-species variance rate
    "haploid effective population size, shared across all populations"
    Ne::Float64
    "λ = v0/(Ne σ²). This is 1 if the root population is at equilibrium"
    lambda::Float64
end
GaussianCoalescent(v0,s2,Ne=1.0) = GaussianCoalescent(v0,s2,Ne,v0/(Ne*bsp_var))
evomodelname(::GaussianCoalescent) = "Gaussian with coalescent"

function gaussiancoalescent_covariancematrix(
    net::HybridNetwork,
    v0::Float64,
    σ2::Float64;
    Ne::Number=1.0,
    checkpreorder::Bool=true,
)
    PN.check_nonmissing_nonnegative_edgelengths(net,
        "The Gaussian with coalescent needs ≥ 0 edge lengths.")
    PN.check_valid_gammas(net,
        "The Gaussian with coalescent needs valid γ inheritances.")
    checkpreorder && preorder!(net)
    # todo: divide edge lengths by Ne to have them in coalescent units
    σ2eq = σ2 * Ne # rate per coal unit = equilibrium within-species variance
    MV = PN.traversal_preorder(net.vec_node,
        init_gaussiancoalmatrix,
        updateroot_gaussiancoalmatrix!,
        updatetree_gaussiancoalmatrix!,
        updatehybrid_gaussiancoalmatrix!,
        σ2eq, v0)
    M = MatrixTopologicalOrder(MV[1], net, :b) # nodes in both columns & rows
    # below: transform V from a vector to a 1-column matrix
    # todo: consider modifying PN.MatrixTopologicalOrder to take an `Array`
    #       more generally, not just matrices?
    V = MatrixTopologicalOrder(reshape(MV[2], (length(net.node),1)), net, :r)
    return M, V
end

function init_gaussiancoalmatrix(nodes::Vector{Node}, params...)
    n = length(nodes)
    M = zeros(Float64,n,n) # (co)variances of species Means
    V = zeros(Float64,n)   # expected within-species Variances
    return([M,V])
end
function updateroot_gaussiancoalmatrix!(MV::Vector, i::Int, bsp_var,v0)
    MV[2][i] = v0
    return true
end

coal_noevent(u)  = exp(-u) # P(no, they do not coalesce)
coal_sharedprop(u) = 1 - (1-exp(u))/u # r(u) = r(l/Ne): proportion shared
coal_sharedtime(u) = u - (1-exp(-u)) # u r(u): time until u shared by coalescence

function updatetree_gaussiancoalmatrix!(
    MV::Vector,
    i::Int,
    parentind::Int,
    edge::Edge,
    σ2, # per coalescent unit: assumes Ne=1 and edge lengths in coalescent units
    args...
)
    M, V = MV
    for j in 1:(i-1)
        M[i,j] = M[j,parentind]
        M[j,i] = M[j,parentind]
    end
    u = edge.length
    p2 = coal_noevent(u)
    p1 = 1-p2
    M[i,i] = σ2 * coal_sharedtime(u) + M[parentind,parentind] + V[parentind] * p1
    V[i] = V[parentind] * p2  +  σ2 * p1
    return true
end
function updatehybrid_gaussiancoalmatrix!(
    MV::Vector,
    i::Int,
    parindx::AbstractVector{Int},
    paredge::AbstractVector{Edge},
    σ2,
    args...
)
    M, V = MV
    for j in 1:(i-1)
        for (pi,pe) in zip(parindx, paredge)
            M[i,j] += pe.gamma * M[pi,j]
        end
        M[j,i] = M[i,j]
    end
    for k1 in eachindex(paredge)
        p1i = parindx[k1]
        p1e = paredge[k1]
        p1u = p1e.length
        p1γ = p1e.gamma
        q1 = 1-coal_noevent(p1u)
        r1 = coal_sharedtime(p1u)
        M[i,i] +=  p1γ^2 * (σ2 * r1 + M[p1i,p1i] + V[p1i] * q1)
        for k2 in (k1+1):length(paredge)
            M[i,i] += 2 * p1γ * paredge[k2].gamma * M[p1i,parindx[k2]]
        end
        V[i] += p1γ * (σ2 * p1u + M[p1i,p1i] + V[p1i])
    end
    V[i] -= M[i,i]
    return true
end
