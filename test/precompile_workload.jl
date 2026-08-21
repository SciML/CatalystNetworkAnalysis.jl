using Catalyst, CatalystNetworkAnalysis, Test

@testset "Precompile workload" begin
    rn = @reaction_network begin
        k₁, A --> B
        k₂, B --> A
    end

    @test isconservative(rn)
    @test isconsistent(rn)
end
