using PhyloTraits
using Test
using Aqua

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
