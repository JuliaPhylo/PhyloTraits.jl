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
    paramlist::Dict,
    xargs...
)
    MV = gaussiancoalescent_covariancematrix_init(net, true, true)
    Ne = gc.Ne
    if Ne != 1
        net = deepcopy(net)
        for e in net.edge  e.length /= Ne; end
    end
    λspec  = paramspec(paramlist, :lambda, lower=0.0)
    v0spec = paramspec(paramlist, :v0, lower=0.0, fixed=false)
    !v0spec.fixed || !λspec.fixed ||
        error("please estimate either λ or v0, for the Gaussian-Coalescent")
    Nespec = paramspec(paramlist, :Ne, lower=0.0, fixed=true)
    Nespec.fixed || error("estimation of Ne is not implemented yet")
    res = phylolm_gcoal_lambda(X,Y,net,reml, MV,λspec,v0spec,
            nonmissing, ind,
            ftolRel, xtolRel, ftolAbs, xtolAbs)
    if Ne != 1
        setNe!(res.evomodel, Ne)
        # sigma2_phylo(res) is still per coalescent units: based on res.lm
    end
    return res
end

function phylolm_gcoal_lambda(
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool,
    MV::Array,
    λspec::ParamSpec,
    v0spec::ParamSpec, # fixit: not used for now. implement method to fix it
    nonmissing::BitArray{1},
    ind::Vector{Int},
    ftolRel::AbstractFloat,
    xtolRel::AbstractFloat,
    ftolAbs::AbstractFloat,
    xtolAbs::AbstractFloat,
)
    (v0spec.fixed && !λspec.fixed) &&
        error("optimization of λ under a fixed v0 has yet to be implemented")
    M, V = gaussiancoalescent_covariancematrix!(MV, net, getvalue(λspec))
    ind_nm = ind[nonmissing] # same length as Y
    Vλ = M[:tips][ind_nm,ind_nm] + Diagonal(V[:tips][ind_nm,1])
    gc_dof = 1     # optimize σ2
    if λspec.fixed # for (each) fixed λ, optimize σ2 = λ*v0 analytically
        λ = getvalue(λspec)
    else
        gc_dof += 1 # also optimize λ
        optsum = OptSummary([getvalue(λspec)],
          [1e-100], # constraint λ≥0 but avoid <0 trials
          :LN_BOBYQA; initial_step=[0.01],
          ftol_rel=ftolRel, ftol_abs=ftolAbs, xtol_rel=xtolRel, xtol_abs=[xtolAbs])
        optsum.maxfeval = 1000 # max number of iterations
        # no upper bound theoretically
        opt = Opt(optsum)
        function fun(x::Vector{Float64}, g::Vector{Float64})
            λ = x[1]
            M, V = gaussiancoalescent_covariancematrix!(MV, net, λ)
            Vλ .= M[:tips][ind_nm,ind_nm] .+ Diagonal(V[:tips][ind_nm,1])
            res = loglik_gcoal_lambda(X,Y,Vλ,reml)
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, optsum.initial)
        λ = xmin[1]
    end
    M, V = gaussiancoalescent_covariancematrix!(MV, net, λ)
    Vλ .= M[:tips][ind_nm,ind_nm] .+ Diagonal(V[:tips][ind_nm,1])
    linmod, Vy, RL, logdetVy = pgls(X,Y,Vλ)
    Ndof = (reml ? dof_residual(linmod) : nobs(linmod))
    σ2eq = deviance(linmod) / Ndof
    v0 = λ * σ2eq
    res = PhyloNetworkLinearModel(linmod, M, Vy, RL, Y, X, logdetVy,
        reml, ind, nonmissing, GaussianCoalescent(v0,σ2eq,1.0,λ, gc_dof))
    return res
end

function loglik_gcoal_lambda(
    X::Matrix,
    Y::Vector,
    Vλ::Matrix,
    reml::Bool;
)
    linmod, Vy, RL, logdetVy = pgls(X,Y,Vλ)
    n = (reml ? dof_residual(linmod) : nobs(linmod))
    res = n*log(deviance(linmod)) + logdetVy
    if reml res += logdet(linmod.pp.chol); end
    return res
end
