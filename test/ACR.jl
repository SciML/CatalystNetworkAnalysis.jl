const C = CatalystNetworkAnalysis

let
    rn = @reaction_network begin
        k1, 2A --> A
        k2, 3A --> 5A
        k3, 4A --> 3A
    end

    p = Dict([k1 => 1., k2 => 1., k3 => 1.])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR

    p = Dict([k1 => 2., k2 => 1.5, k3 => 1.])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR

    # This reaction network has four ACR regions. 
    rn = @reaction_network begin
        @parameters α β
        (2α + β), A --> 0
        (3 + αβ, α + β), A <--> 2A
        3, 2A --> 3A
        1, 3A --> 4A
    end
    
    p1 = Dict([α => 1.5, β => 2.6])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR
    # this case has no acr
    p2 = Dict([α => 1.5, β => 2.5])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR
    p3 = Dict([α => 1.5, β => 1.6])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR
    p4 = Dict([α => 0.5, β => 2.6])
    @test C.isconcentrationrobust(rn, p) == :MASS_ACTION_ACR


    # This reaction network has unconditional ACR, but is not detected
    rn = @reaction_network begin
        k1, B --> A
        k2, 2A + B --> A + 2B
    end
    @test C.isconcentrationrobust(rn) == :INCONCLUSIVE
end

# Examples with detection
let
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
end
