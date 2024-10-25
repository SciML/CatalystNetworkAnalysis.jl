const C = CatalystNetworkAnalysis

# Basic examples (section 2) 
let
    # Detected by the deficiency one check.
    rn = @reaction_network begin
        k1, B --> A
        k2, 2A + B --> A + 2B
    end

    p = Dict([:k1 => 1, :k2 => 1])
    @test C.isconcentrationrobust(rn; p = p) == :GLOBAL_ACR

    # This reaction network has four ACR regions. 
    rn = @reaction_network begin
        @parameters α β
        (2α + β), A --> 0
        (3 + α*β, α + β), A <--> 2A
        3, 2A --> 3A
        1, 3A --> 4A
    end
    
    p1 = Dict([:α => 1.5, :β => 2.6])
    @test C.isconcentrationrobust(rn; p = p1) == :INCONCLUSIVE
    # this case has no acr
    p2 = Dict([:α => 1.5, :β => 2.5])
    @test C.isconcentrationrobust(rn; p = p2) == :INCONCLUSIVE
    p3 = Dict([:α => 1.5, :β => 1.6])
    @test C.isconcentrationrobust(rn; p = p3) == :INCONCLUSIVE
    p4 = Dict([:α => 0.5, :β => 2.6])
    @test C.isconcentrationrobust(rn; p = p4) == :INCONCLUSIVE
end

# Examples with detection
let
    # Via Saturation
    rn = @reaction_network begin
        (k1, k2), S1 + E <--> C1
        k3, C1 --> S2 + E
        (k4, k5), S2 + E <--> C2
        k6, C2 --> S3 + E
        (k7, k8), S2 + C1 <--> C3
        k9, C3 --> S3 + C1
        (k10, k11), S3 + C1 <--> C4
        k12, C4 --> S1 + C1
    end
    @test C.isconcentrationrobust(rn) == :GLOBAL_ACR

    # Via elimination ideals
    rn = @reaction_network begin
        k1, 0 --> A
        k2, 2A --> 0
        k3, 2B --> 3B
        k4, A + B --> A
        k5, B --> 2B
    end
    p = Dict([:k1 => 2, :k2 => 1, :k3 => 1, :k4 => 1, :k5 => 1.41])
    @test C.isconcentrationrobust(rn; p = p) == :MASS_ACTION_ACR
end

# Testing the utilities for the necessary condition on absolute concentration robustnesss.
let
    rn = @reaction_network begin
        (k1, k2), 2A <--> B
        k3, B --> C
        (k4, k5), B + C <--> D
        k6, D --> 2B
        k7, 2B --> A + E
        k8, A + E --> F
        (k9, k10), 2B <--> F
    end
    @test deficiency(rn) == 1

    S, D = CatalystNetworkAnalysis.removespec(rn, 2)
    @test S == [-2 2 0 0 0 0 1 -1 0 0;
                 0 0 1 -1 1 0 0 0 0 0; 
                 0 0 0 1 -1 -1 0 0 0 0;
                 0 0 0 0 0 0 1 -1 0 0;
                 0 0 0 0 0 0 0 1 1 -1]

    @test CatalystNetworkAnalysis.deficiency(S, D) == 0

    S, D = CatalystNetworkAnalysis.removespec(rn, 3)
    @test S == [-2 2 0 0 0 0 1 -1 0 0;
                 1 -1 -1 -1 1 2 -2 0 -2 2;
                 0 0 0 1 -1 -1 0 0 0 0;
                 0 0 0 0 0 0 1 -1 0 0;
                 0 0 0 0 0 0 0 1 1 -1]

    @test CatalystNetworkAnalysis.deficiency(S, D) == 1
end
