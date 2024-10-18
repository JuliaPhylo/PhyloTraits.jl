using PhyloTraits
using Test
using Aqua

using BioSymbols
using CSV
using DataFrames
using GLM # for coef, nobs, residuals etc.
using LinearAlgebra: norm, diag, logdet, PosDefException # LinearAlgebra.rotate! not brought into scope
using PhyloNetworks
using Random
using StableRNGs
using StaticArrays # for rate substitution matrices
using Statistics
using StatsBase # for aic etc., stderr
using StatsAPI

@testset "PhyloTraits Code quality (Aqua.jl)" begin
    Test.detect_ambiguities(PhyloTraits)
    Aqua.test_all(
        PhyloTraits;
        ambiguities = (broken=false),
        persistent_tasks = false,
    )
end
@testset "PhyloTraits.jl" begin
    include("test_lm.jl")
    include("test_lm_tree.jl")
    include("test_lm_withinspecies.jl")
    include("test_traits_discrete.jl")
    include("test_simulate.jl")
    include("test_simulate_mbd.jl")
end
