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
    if gc.Ne != 1
        net = deepcopy(net)
        for e in net.edge  e.length /= Ne; end
    end
    v0spec = paramspec(paramlist, :v0, lower=0.0)
    !v0spec.fixed ||
        error("for now, the Gaussian-Coalescent model cannot fix v0")
    λspec = paramspec(paramlist, :lambda, lower=0.0)
    # fixit: synchronize σ2eqm, v0spec and λspec later, when several are used
    phylolm_gcoal_lambda(X,Y,net,reml, MV,λspec,v0spec,
            nonmissing, ind,
            ftolRel, xtolRel, ftolAbs, xtolAbs,
            fix_v0)
end

function phylolm_gcoal_lambda(
    X::Matrix,
    Y::Vector,
    net::HybridNetwork,
    reml::Bool,
    MV::Array,
    λspec::ParamSpec,
    v0spec::Float64,
    nonmissing::BitArray{1},
    ind::Vector{Int},
    ftolRel::AbstractFloat,
    xtolRel::AbstractFloat,
    ftolAbs::AbstractFloat,
    xtolAbs::AbstractFloat,
)
    M, V = gaussiancoalescent_covariancematrix!(MV, net, getvalue(λspec))
    ind_nm = ind[nonmissing] # same length as Y
    Vind = M[:tips][ind_nm,ind_nm]
    @show Vind
    Vind += Diagonal(V[:tips][ind_nm,1])
    @show Vind
    if λspec.fixed
        λ = getvalue(λspec)
        @warn "fixed lambda=$λ for now"
        linmod, Vy, RL, logdetVy = pgls(X,Y,Vind)
        res = PhyloNetworkLinearModel(linmod, M, Vy, RL, Y, X, logdetVy,
            reml, ind, nonmissing, GaussianCoalescent(v0, σ2eq))
        return res

        opt = NLopt.Opt(:LN_BOBYQA, 1)
        NLopt.ftol_rel!(opt, ftolRel) # relative, objective
        NLopt.ftol_abs!(opt, ftolAbs) # absolute, objective
        NLopt.xtol_rel!(opt, xtolRel) # relative, parameter
        NLopt.xtol_abs!(opt, xtolAbs) # absolute, parameter
        NLopt.maxeval!(opt, 1000) # max number of iterations
        NLopt.lower_bounds!(opt, 1e-100) # constraint λ≥0 but avoid <0 trials
        # no upper bound theoretically
        function fun(x::Vector{Float64}, g::Vector{Float64})
            x = convert(AbstractFloat, x[1])
            res = loglik_gcoal_lambda(x, X,Y,V, reml, gammas, times; nonmissing=nonmissing, ind=ind)
            # count =+ 1
            #println("f_$count: $(round(res, digits=5)), x: $(x)")
            return res
        end
        NLopt.min_objective!(opt, fun)
        fmin, xmin, ret = NLopt.optimize(opt, [startingValue])
        # Best value dans result
        res_lam = xmin[1]

    end
    error("optimization of v0 has yet to be implemented")
end

function loglik_gcoal_lambda(
    lam::AbstractFloat,
    X::Matrix,
    Y::Vector,
    Vind::Matrix,
    reml::Bool;
)
    reml && error("ML implemented only so far. Use reml=false")
    linmod, Vy, RL, logdetVy = pgls(X,Y,Vind)
    n = nobs(linmod)
    # n = (reml ? dof_residual(linmod) : nobs(linmod))
    res = n*log(deviance(linmod)) + logdetVy
    # if reml res += logdet(linmod.pp.chol); end
    return res
end
