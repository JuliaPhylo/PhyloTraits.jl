function phylolm(
    gc::GaussianCoalescent,
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
    startingValue::Real, # σ2eq, not used: taken from gc
    fixedValue::Union{Real,Missing}=missing, # v0
    xargs...
)
    MV = gaussiancoalescent_covariancematrix_init(net, true, true)
    if gc.Ne != 1
        net = deepcopy(net)
        for e in net.edge  e.length /= Ne; end
    end
    σ2eq = gc.sigma2 * gc.Ne # equilibrium within-species variance
    fix_v0 = !ismissing(fixedValue)
    fix_v0 || error("for now, the Gaussian-Coalescent model requires a fixed v0")
    v0 = fixedValue
    phylolm_gcoal_sigma2(X,Y,net,reml, MV,v0,σ2eq,
            nonmissing, ind,
            ftolRel, xtolRel, ftolAbs, xtolAbs,
            fix_v0)
end

function phylolm_gcoal_sigma2(
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool,
    MV::Array,
    v0::Float64,
    σ2eq::Float64,
    nonmissing::BitArray{1},
    ind::Vector{Int},
    ftolRel::AbstractFloat,
    xtolRel::AbstractFloat,
    ftolAbs::AbstractFloat,
    xtolAbs::AbstractFloat,
    fix_v0::Bool,
)
    if fix_v0
        @warn "fixed sigma2=$σ2eq for now"
        M, V = gaussiancoalescent_covariancematrix!(MV, net, v0, σ2eq)
        ind_nm = ind[nonmissing] # same length as Y
        Vind = M[:tips][ind_nm,ind_nm]
        @show Vind
        Vind += Diagonal(V[:tips][ind_nm,1])
        @show Vind
        linmod, Vy, RL, logdetVy = pgls(X,Y,Vind)
        res = PhyloNetworkLinearModel(linmod, M, Vy, RL, Y, X, logdetVy,
            reml, ind, nonmissing, GaussianCoalescent(v0, σ2eq))
        return res
    end
    error("optimization of v0 has yet to be implemented")
end
