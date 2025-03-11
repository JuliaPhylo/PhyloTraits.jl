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
    predict(
        obj::ReconstructedStates;
        interval::Union{Symbol,Nothing}=nothing,
        level::Real=0.95,
        text::Bool=false,
        digits::Int=2,
        missingmark::AbstractString="*",
        combine::Bool=false
    )

Estimated reconstructed states and, if `interval=:prediction`, prediction
intervals with level `level` at all internal nodes and tips with missing traits.

If `text=true`, the prediction and intervals are formated as string for easy
plotting with the `nodelabel` argument to `plot` from
package [`PhyloPlots`](https://github.com/juliaphylo/PhyloPlots.jl).
In that case,
`digits` controls the number of digits shown,
`missingmark` adds a distinctive mark to the prediction of tips with missing data,
(set to `missingmark=""` for no mark),
and if `combine=true`, the prediction and bound of the intervals are combined
into a single text string.
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
All methods use the `showWarnings` argument to specify whether warnings are printed
core method called by every other method:

1. ancestralreconstruction(Vzz, VzyVyinvchol, RL, Y, m_y, m_z,
                                nodenumbers, tipnumbers, sigma2, add_var, model)

higher-level methods, for real data:

2. ancestralreconstruction(dataframe, net; tipnames=:tipnames, showWarnings, kwargs...)
   - dataframe: 2 columns only, species names & tip response values
   - fits an intercept-only model, then calls #3
   - by default without kwargs: model = BM w/o within-species variation

3. ancestralreconstruction(PhyloNetworkLinearModel[, Matrix]; showWarnings)
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
    model::Union{PhyloNetworkLinearModel,Missing}=missing)
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
function ancestralreconstruction(obj::PhyloNetworkLinearModel, X_n::Matrix; showWarnings::Bool=true)
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
        showWarnings && @warn """There were no indication for the position of the tips on the network.
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
    showWarnings && @warn """These prediction intervals show uncertainty in ancestral values,
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
    ancestralreconstruction(obj::PhyloNetworkLinearModel[, X_n::Matrix]; showWarnings::Bool=true)

Function to find the ancestral traits reconstruction on a network, given an
object fitted by function [`phylolm`](@ref). By default, the function assumes
that the regressor is just an intercept. If the value of the regressor for
all the ancestral states is known, it can be entered in X_n, a matrix with as
many columns as the number of predictors used, and as many lines as the number
of unknown nodes or tips.
Setting the optional keyword argument `showWarnings` to `false` supresses
PhyloTraits-specific warning outputs.

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
function ancestralreconstruction(obj::PhyloNetworkLinearModel; showWarnings::Bool=true)
    # default reconstruction for known predictors
    if ((size(obj.X)[2] != 1) || !any(obj.X .== 1)) # Test if the regressor is just an intercept.
        error("""Predictor(s) other than a plain intercept are used in this `PhyloNetworkLinearModel` object.
    These predictors are unobserved at ancestral nodes, so they cannot be used
    for the ancestral state reconstruction. If these ancestral predictor values
    are known, please provide them as a matrix argument to the function.
    Otherwise, you might consider doing a multivariate linear regression (not implemented yet).""")
    end
    X_n = ones((length(obj.V.nodenumbers_toporder) - sum(obj.nonmissing), 1))
    ancestralreconstruction(obj, X_n; showWarnings)
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
    showWarnings::Bool=true,
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
    return ancestralreconstruction(reg; showWarnings)
end
