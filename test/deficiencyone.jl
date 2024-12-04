# Testing whether regularity is properly tested. 

# This file contains tests for the implementation of the deficiency one algorithm and the higher deficiency algorithm. 
using Catalyst
import CatalystNetworkAnalysis as C
using JuMP, HiGHS

let
    regular_rn = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        (k8, k9), C <--> D
        k10, D --> A
        k11, D --> B + E
        k12, B + E --> C
    end

    # Fails the terminal linkage class cut-link condition. 
    irregular_rn = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        (k8, k9), C <--> D
        k10, D --> A
        k11, D --> B + E
        k12, B + E --> C
        k13, E --> F
    end

    # Fails the number of terminal linkage class condition. D --> C is now irreversible, and C <--> B + E now is, making (C, B + E) a terminal linkage class. 
    irregular_rn2 = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        k8, D --> C
        k9, D --> A
        k10, D --> B + E
        (k11, k12), B + E <--> C
    end

    reactioncomplexes(regular_rn)
    reactioncomplexes(irregular_rn)
    reactioncomplexes(irregular_rn2)
    @test C.isregular(regular_rn) == true
    @test C.isregular(irregular_rn) == false
    @test C.isregular(irregular_rn2) == false
end

# Testing `generatepartitions()`: whether the number of partitions matches expectations

let
    # This reaction network has {A} as a trivial terminal linkage class. 
    one_trivial_tlc = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        (k8, k9), C <--> D
        k10, D --> A
        k11, D --> B + E
        k12, B + E --> C
    end

    # By making D <--> A reversible, the terminal linkage class containing A gets expanded. 
    zero_trivial_tlc = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        (k8, k9), C <--> D
        (k10, k11), D <--> A
        k12, D --> B + E
        k13, B + E --> C
    end

    # By making D --> C irreversible, we add another terminal linkage class, {C}
    multiple_trivial_tlcs = @reaction_network begin
        k1, 2A --> B
        (k2, k3), B <--> C + D
        (k4, k5), C + D <--> F
        (k6, k7), C + D <--> E
        k8, D --> C
        k9, D --> A
        k10, D --> B + E
        k11, B + E --> C
    end

    part1 = C.generatepartitions(one_trivial_tlc)
    part0 = C.generatepartitions(zero_trivial_tlc)
    partn = C.generatepartitions(multiple_trivial_tlcs)

    @test length(part1) == 1 # (3^1+1) / 2 - 1
    @test length(part0) == 1 # (3^2+1) / 2 - 2^2
    @test length(partn) == 2 # (3^1+1) / 2
end

# Testing whether the proper confluence vector orientation is identified
# and generates a solution

let
    rn = @reaction_network begin
        (k1, k2), A + S <--> AS
        (k3, k4), B + S <--> BS
        k5, AS + BS --> C + 2S
        (k6, k7), A <--> 0
        (k8, k9), B <--> 0
        k10, C --> 0
    end

    part = C.generatepartitions(rn)
    @test length(part) == 13 # (3^3 + 1)/2 - 1
    g = C.confluencevector(rn)
    @test C.confluencevector(rn) ≈ [-1., 1., -1., 1., -1., 1., 1., -1., 1., -1.]
    # U, M, L
    
    correctpartition = [[1, 2, 7, 8, 9], [5, 10], [3, 4]]
    cutdict = C.cutlinkpartitions(rn)

    Y = complexstoichmat(rn); S = netstoichmat(rn)
    s, c = size(Y); r = size(S, 2)
    
    # Initialization
    model = Model(HiGHS.Optimizer); set_silent(model)
    @variable(model, μ[1:s]); @variable(model, n)
    @objective(model, Min, 0)

    μ_sol, feasible = C.solveconstraints(rn, model, g, correctpartition, cutdict)
    @test feasible == true
end

# Testing whether the deficiency one algorithm returns the correct answer
# in some simple cases. 

let
    rn1 = @reaction_network begin
        (k1, k2), A + B <--> 2A
        (k3, k4), A <--> 0
        (k5, k6), B <--> 0
    end

    rn2 = @reaction_network begin
        (k1, k2), 2A + B <--> 3A
        (k3, k4), A <--> 0
        (k5, k6), B <--> 0
    end

    rn3 = @reaction_network begin
        (k1, k2), 2A + B <--> 3A
        (k3, k4), A <--> 0
        (k5, k6), B <--> A
    end

    # Lotka-Volterra network. Has periodic solutions, but only one equilibrium
    # per stoichiometric compatibility class. 
    rn4 = @reaction_network begin
        k1, A --> 2A
        k2, A + B --> 2B
        k3, B --> 0
    end

    # Networks that can admit multiple equilibria. 
    # Edelstein network. 
    rn5 = @reaction_network begin
        (k1, k2), A <--> 2A
        (k3, k4), A + B <--> C
        (k5, k6), B <--> C
    end

    # Chirality pattern formation network 
    rn6 = @reaction_network begin
        (k1, k2), L + 2R + P <--> 3R + Q
        (k3, k4), R + 2L + P <--> 3L + Q
        (k5, k6), P <--> 0
        (k7, k8), 0 <--> Q
    end
    
    @test all(C.isregular, [rn1, rn2, rn3]) == true
    @time @test C.deficiencyonealgorithm(rn1) == false
    @time @test C.deficiencyonealgorithm(rn2) == true 
    @time @test C.deficiencyonealgorithm(rn3) == false
    @time @test C.deficiencyonealgorithm(rn4) == false 
    @time @test C.deficiencyonealgorithm(rn5) == true 
    @time @test C.deficiencyonealgorithm(rn6) == true 
end


### Higher Deficiency Algorithm
let
    rn = @reaction_network begin
        (k1, k2), E1 + S1 <--> E1S1
        (k3, k4), E1S1 <--> E1 + S2 
        (k5, k6), E1 + S2 <--> E1S2
        k7, E1S2 --> E1 + S3

        (g1, g2), E2 + S3 <--> E2S3
        g3, E2S3 --> E2 + S2 
        (g5, g6), E2 + S2 <--> E2S2
        g7, E2S2 --> E2 + S1
        
        (f1, f2), E3 + S1 <--> E3S1
        f3, E3S1 --> E3 + S3 
        (f5, f6), E3 + S3 <--> E3S3
        f7, E3S3 --> E3 + S2
    end
end
