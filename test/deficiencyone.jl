# Testing whether regularity is properly tested. 

import CatalystNetworkAnalysis as C

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

# Testing `generatepartitions()`: whether the number of partitions matches 
# expectations

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

    @test length(part1) == 4 # (3^2+1) / 2 - 1
    @test length(part0) == 1 # (3^2+1) / 2 - 2^2
    @test length(partn) == 5 # (3^3+1) / 2
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

    C.confluencevector(rn)
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

    @test all(C.isregular, [rn1, rn2, rn3]) == true
    @test C.deficiencyonealgorithm(rn1) == false
    @test C.deficiencyonealgorithm(rn2) == true 
    @test C.deficiencyonealgorithm(rn3) == false
end
