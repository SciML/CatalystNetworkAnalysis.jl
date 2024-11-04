include("../test_networks.jl")
const C = CatalystNetworkAnalysis

# Tests for the NetworkSummary function.

# Ground truths
networks = [MAPK, zigzag, EnvZ_OmpR,
            his_kinase, oneSitePhosphorylation(10),
            cellDeathNetwork(10), edelstein(10)]

# Testing bistability is detected
let
    for (name, network) in nets
        @unpack rn, u0, p = network
        ns = C.networksummary(rn; u0 = u0, p = p)
        @test ns.steadystates == :INCONCLUSIVE || :KINETICALLY_MULTIPLE || :GLOBALLY_MULTIPLE
    end
end

