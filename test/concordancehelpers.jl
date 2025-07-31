import CatalystNetworkAnalysis as C
using Catalyst
using StableRNGs

rng = StableRNG(123)
# moveforward and movebackward
let
    sp = [1, 0, 1, 0, 1, 0]
    fixedsigns = [6, 7, 8]
    n = 8
    @test C.moveforward(sp, fixedsigns, n) == [1, 0, 1, 0, 1, 0, 0, 0]
    @test C.movebackward(sp, fixedsigns) == [1, 0, 1, 0, -1]
    @test C.movebackward(sp, fixedsigns) == [1, 0, 1, 0, 0]

    sp2 = Int[]
    @test C.moveforward(sp2, fixedsigns, n) == [1]
    @test C.movebackward(sp2, fixedsigns) == [0]
    @test C.movebackward(sp2, fixedsigns) == []
end

# sign compatibility tests
let
    S = [1 1 0;
         0 1 0;
         1 0 0]
    n = 3

    @test C.issigncompatible(S, Float64[], freeindices = collect(1:n)) == true
    @test_throws "The number of free signs and assigned signs does not sum to the length of the vector." C.issigncompatible(
        S, Float64[])
    @test C.issigncompatible(S, [1, 1, 1]) == true
    @test C.issigncompatible(S, [-1, 1, 1]) == false
end

# generating α sign patterns
# Concordance tests
let
    rn = @reaction_network begin
        (k1, k2), A + B <--> P
        (k3, k4), B <--> 2A
    end

    rn2 = @reaction_network begin
        (k1, k2), A + B <--> P
        (k3, k4), B + C <--> Q
        (k5, k6), C <--> 2A
    end

    rn3 = @reaction_network begin
        (k1, k2), A + B <--> P
        (k3, k4), B + C <--> Q
        (k5, k6), C + D <--> R
        (k7, k8), D <--> 2A
    end

    σ_1 = [1, 0, -1]
    #  A, B, P, C, Q
    σ_2 = [0, 1, 1, 0, -1]
    #  A, B, P, C, Q,  D, R
    σ_3 = [0, 0, 1, 1, -1, -1, 1]

    @test C.generate_α_signpattern(rn, σ_1) == ([1, -1, 0, 1], Int[])
    @test C.generate_α_signpattern(rn2, σ_2) == ([1, 1, 1, -1, 0, 0], Int[])
    @test C.generate_α_signpattern(rn3, σ_3) == ([0, 1, 1, -1, 1, -1, 0], [5])

    @test C.isconcordant(rn) == true
    @test C.isconcordant(rn2) == false
    @test C.isconcordant(rn3) == true
end

let
    WNT = @reaction_network begin
        k1, A10 --> ∅
        (k2, k3), ∅ <--> A11
        (k4, k5), ∅ <--> A12
        (k6, k7), A1 <--> A2
        k8, A2 + A4 --> A16
        k9, A16 --> A2 + A6
        (k10, k11), A3 <--> A4
        (k12, k13), A4 <--> A6
        (k14, k15), A6 <--> A7 + A12
        (k16, k17), A3 + A11 <--> A8
        k18, A8 --> A9
        k19, A9 --> A3 + A10
        (k20, k21), A11 + A13 <--> A14
    end

    @test C.isconcordant(WNT) == true
end
