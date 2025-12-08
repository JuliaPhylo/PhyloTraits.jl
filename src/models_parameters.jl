struct ParamSpec{T<:AbstractFloat}
    "value to start optimization, or to get fixed"
    start::T
    "lower bound"
    lower::T
    "upper bound"
    upper::T
    "should this parameter be fixed (not optimized)?"
    fixed::Bool
end
function ParamSpec(
    x::NamedTuple,
    d::ParamSpec{T}=ParamSpec{Float64}(1.0, -Inf, Inf, false)
) where T
    s = get(x, :start, missing)
    l = get(x, :lower, missing)
    u = get(x, :upper, missing)
    f = get(x, :fixed, missing)
    if !ismissing(l) && !ismissing(u) && l==u
        f = true
    end
    if ismissing(l) l = d.lower; end
    if ismissing(u) u = d.upper; end
    l ≤ u || error("lower bound $l should be ≤ upper bound $u")
    if ismissing(s)
        s = max(d.start, l) # start = lower if 'start' not provided but 'lower' is
    end
    if ismissing(f) f = d.fixed; end
    f || # case where the parameter is not fixed
        (l ≤ s ≤ u) || # check inclusion
        error("starting value $s needs to be in between bounds [$l,$u]")
    !f || # case where the parameter is fixed
        (l ≤ s ≤ u) || # check inclusion
        (s == 0.0 && isapprox(l, s; atol = eps(typeof(l)))) || # case where lower bound is 1e-100
        error("fixed value $s needs to be in between bounds [$l,$u]")
    return ParamSpec{T}(s,l,u,f)
end

getvalue(p::ParamSpec) = p.start

function paramspec(
    paramlist::Dict,
    paramname::Symbol;
    start=1.0,
    lower=-Inf,
    upper=Inf,
    fixed=false,
)
    default = ParamSpec{Float64}(start, lower, upper, fixed)
    if haskey(paramlist, paramname)
        return ParamSpec(paramlist[paramname], default)
    end
    return default
end
