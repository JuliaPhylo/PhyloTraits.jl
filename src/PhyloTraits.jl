module PhyloTraits

# default tolerances to optimize parameters in continuous trait evolution models
# like lambda, sigma2_withinspecies / sigma2_BM, etc.
const fAbsTr = 1e-10
const fRelTr = 1e-12
const xAbsTr = 1e-10
const xRelTr = 1e-10
const alphaRASmin = 0.02
const alphaRASmax = 50.0
const pinvRASmin = 1e-8
const pinvRASmax = 0.99
const kappamax = 20.0

using BioSequences
using BioSymbols
using DataFrames # innerjoin new in v0.21
using Distributions #for RateVariationAcrossSites
using FASTX
using GLM
using LinearAlgebra: diag, I, logdet, norm, LowerTriangular, mul!, lmul!, rmul!,
        Diagonal, cholesky, qr, BLAS
# caution: both LinearAlgebra and PhyloNetworks export rotate!
# alternative: drop support for julia v1.4, as LinearAlgebra.rotate! requires julia v1.5
# using LinearAlgebra # bring all of LinearAlgebra into scope
# import LinearAlgebra.rotate! # allow re-definition of rotate!
using NLopt
using PhyloNetworks
using PhyloNetworks: Edge, Node, MatrixTopologicalOrder
using Printf: @printf, @sprintf
using Random
using Random: AbstractRNG, default_rng
using StaticArrays
using Statistics: mean, quantile, median
using StatsAPI: StatsAPI, coef, coefnames, coeftable, confint, deviance
using StatsAPI: dof, dof_residual, fit, fit!, fitted, isfitted, islinear, leverage
using StatsAPI: loglikelihood, modelmatrix, nobs, predict, r2, residuals
using StatsAPI: response, stderror, vcov
using StatsBase
using StatsFuns # logsumexp, logaddexp, log2π, various cdf
using StatsModels # re-exported by GLM. for ModelFrame ModelMatrix Formula etc

const PN = PhyloNetworks

# import: to extend methods from othe packages with new methods defined here
import Base: show
import GLM: ftest, fit!
import PhyloNetworks: tiplabels
import StatsModels: coefnames

export ftest # from GLM
# continuous traits
export phylolm, PhyloNetworkLinearModel
export simulate, TraitSimulation
export ParamsBM, ParamsMultiBM
export ShiftNet, shiftHybrid, getShiftEdgeNumber, getShiftValue
export regressorShift, regressorHybrid
export ancestralStateReconstruction, ReconstructedStates
export sigma2_phylo, sigma2_within
export mu_phylo
export lambda_estim
export expectations, expectationsPlot
export predint, predintPlot
# discrete traits
export TraitSubstitutionModel
export EqualRatesSubstitutionModel, BinaryTraitSubstitutionModel
export TwoBinaryTraitSubstitutionModel
export JC69, HKY85
export nstates
export Q
export getlabels
export nparams
export RateVariationAcrossSites
export randomTrait, randomTrait!
export fitdiscrete
export readfastatodna
export stationary
export empiricalDNAfrequencies

include("nloptsummary.jl")
include("traits_continuous.jl")
include("simulate_continuous.jl")
include("substitutionmodels.jl")
include("traits_discrete.jl")

end
