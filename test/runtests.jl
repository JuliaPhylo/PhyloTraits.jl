using PhyloTraits
using Test
using Aqua

@testset "PhyloTraits.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(PhyloTraits)
    end
    # Write your tests here.
end
