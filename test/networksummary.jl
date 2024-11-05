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


# Run on the SBML Mass Action Networks
let
    for (name, network) in ma_SBMLnets
        @unpack rn, u0, p = network
        ma = all(r -> ismassaction(r, rn), reactions(rn))
        @btime ns = C.networksummary(rn; u0 = u0, p = p)
        show(ns)
    end
end 

