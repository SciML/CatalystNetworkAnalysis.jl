using CatalystNetworkAnalysis, Catalyst

# Check if persistence(rn) correctly identifies the following networks as persistent
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

    @test ispersistent(rn) == true
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
    
    
    @test ispersistent(rn) == true
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

    @test ispersistent(rn) == true
end

# The following example cannot be determined to be persistent from the persistence algorithm.

let
    rn = @reaction_network begin
        k1, 2A + B --> C
        k2, C --> A + 2B
        k3, A + 2B --> D
        k4, D --> 2A + B
    end

    @test_throws ErrorException ispersistent(rn)
end
