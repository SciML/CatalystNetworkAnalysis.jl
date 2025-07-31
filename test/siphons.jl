using CatalystNetworkAnalysis, Catalyst

# Check if minimalsiphons() correctly identifies the set of siphons for each of the following reaction networks. 
let
    rn = @reaction_network begin
        (k1, k2), 2A + C <--> A + D
        (k3, k4), A + D <--> E
        (k5, k6), E <--> B + C
        (k7, k8), B + C <--> 2A + C
    end

    siphons = [[1, 4, 5], [1, 2, 4], [2, 3, 4]]
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        (k1, k2), S + E <--> Q
        (k3, k4), Q <--> P + E
        (k5, k6), Q + I <--> R
    end

    siphons = [[2, 3, 6], [5, 6], [1, 3, 4, 6]]
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        (k1, k2), S0 + E <--> X
        k3, X --> P + E
        (k5, k6), P + F <--> Y
        k7, Y --> S0 + F
    end

    siphons = [[2, 3], [5, 6], [1, 3, 4, 6]]
    S = netstoichmat(rn)
    @test all(s->!iscritical(s, S), siphons) == true
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        (k1, k2), E + S0 <--> ES0
        k3, ES0 --> E + S1
        (k4, k5), E + S1 <--> ES1
        k6, ES1 --> E + S2
        (k7, k8), F + S2 <--> FS2
        k9, FS2 --> F + S1
        (k10, k11), F + S1 <--> FS1
        k12, FS1 --> F + S0
    end

    siphons = [[1, 3, 5], [7, 8, 9],
        [2, 3, 4, 5, 6, 8, 9]]
    S = netstoichmat(rn)
    @test all(s->!iscritical(s, S), siphons) == true
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        (k1, k2), E + S0 <--> ES0
        k3, ES0 --> E + S1
        (k4, k5), E + S1 <--> ES1
        k6, ES1 --> E + E_
        (k7, k8), F + E_ <--> FS2
        k9, FS2 --> F + S1
        (k10, k11), F + S1 <--> FS1
        k12, FS1 --> F + S0
        (k13, k14), E_ + S0_ <--> ES0_
        k15, ES0_ --> E_ + S1_
        (k16, k17), E_ + S1_ <--> ES1_
        k18, ES1_ --> E_ + S2_
        (k19, k20), F_ + S2_ <--> FS2_
        k21, FS2_ --> F_ + S1_
        (k22, k23), F_ + S1_ <--> FS1_
        k24, FS1_ --> F_ + S0_
    end

    siphons = [[1, 3, 5],
        [7, 8, 9],
        [15, 16, 17],
        [10, 11, 12, 13, 14, 16, 17],
        [2, 3, 4, 5, 6, 8, 9, 11, 13]]

    S = netstoichmat(rn)
    @test all(s->!iscritical(s, S), siphons) == true
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        (k1, k2), M + E <--> ME
        k3, ME --> My + E
        (k4, k5), My + E <--> MyE
        k6, MyE --> M2 + E
        (k7, k8), M + E <--> ME_
        k9, ME_ --> Mt + E
        (k10, k11), Mt + E <--> MtE
        k12, MtE --> M2 + E
        (k13, k14), M2 + F <--> M2F
        k15, M2F --> My + F
        (k16, k17), My + F <--> MyF
        k18, MyF --> M + F
        (k19, k20), M2 + F <--> M2F_
        k21, M2F_ --> Mt + F
        (k22, k23), Mt + F <--> MtF
        k24, MtF --> M + F
    end

    siphons = [[2, 3, 5, 7, 9],
        [10, 11, 12, 13, 14],
        [1, 3, 4, 5, 6, 7, 8, 9, 11, 12, 13, 14]]
    S = netstoichmat(rn)
    @test all(s->!iscritical(s, S), siphons) == true
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end

let
    rn = @reaction_network begin
        k1, 2A + B --> C
        k2, C --> A + 2B
        k3, A + 2B --> D
        k4, D --> 2A + B
    end

    siphons = [[1, 3, 4], [2, 3, 4]]
    S = netstoichmat(rn)

    # Both of the siphons in this system are critical. In the previous examples, all of the siphons were not critical.  
    @test all(s->iscritical(s, S), siphons) == true
    @time ms = minimalsiphons(rn)
    @test issetequal(ms, siphons)
end
