include(joinpath(@__DIR__, "test_networks.jl"))
using BenchmarkTools
using StableRNGs
const C = CatalystNetworkAnalysis
import HomotopyContinuation

rng = StableRNG(444)

# Correctness tests for the NetworkSummary function for a large set of mass action SBML networks
# Ground truths
networks = [MAPK, zigzag, EnvZ_OmpR,
    his_kinase, oneSitePhosphorylation(10),
    cellDeathNetwork(10), edelstein(10)]

# Test the correctness of each function.
# 1) For concentration-robust networks, find a witness for concentration robustness
# 2) For steady states, do steady state computations
# 3) For persistence
num_specs = Int64[]
num_rxs = Int64[]
ns_times = Float64[]
hc_times = Float64[]

num_acr = 0
num_mss = 0
num_acr_detected = 0
num_mss_detected = 0

for (name, network) in ma_nets
    @unpack rn, u0, p = network
    println("Computing network summary for $name")
    @btime ns = C.networksummary(rn; u0, p)
    push!(ns_times, ns)
    push!(num_rxs, length(reactions(rn)))
    push!(num_specs, length(species(rn)))

    println("Computing steady states for $name")
    nss = hc_steady_states(rn, ps = p; u0)
    if ns.steadystates == :STRUCTURALLY_UNIQUE || :KINETICALLY_UNIQUE
        @test length(nss) == 1
    elseif ns.steadystates == :STRUCTURALLY_MULTIPLE || :KINETICALLY_MULTIPLE
        @test length(nss) > 1
        length(nss) > 1 && begin
            num_mss_detected += 1
            num_mss += 1
        end
    else
        length(nss) > 1 && (num_mss += 1)
    end

    if ns.concentrationrobust == :GLOBAL_ACR || :MASS_ACTION_ACR
        rs = rn.robustspecies
        nl1 = NonlinearProblem(rn, u_guess = u0, p)
        ss1 = sol(nl1)

        randinit = rand(length(u0))
        randu0 = Dict(zip(species(rn), randinit))
        nl2 = NonlinearProblem(rn, u_guess = randu0, p)
        ss2 = sol(nl2)

        @test ss1 != ss2
        @test ss1.u[rs] ≈ ss2[rs]
        num_acr_detected += 1
    end
end
