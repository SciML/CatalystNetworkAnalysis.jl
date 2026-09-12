using BenchmarkTools
using Catalyst
using CatalystNetworkAnalysis
const C = CatalystNetworkAnalysis

const SUITE = BenchmarkGroup()

rn_simple = @reaction_network begin
    k1, A --> B
    k2, B --> A
end

rn_michaelis_menten = @reaction_network begin
    (k1, k2), S + E <--> C
    k3, C --> P + E
end

rn_acr = @reaction_network begin
    (k1, k2), S1 + E <--> C1
    k3, C1 --> S2 + E
    (k4, k5), S2 + E <--> C2
    k6, C2 --> S3 + E
    (k7, k8), S2 + C1 <--> C3
    k9, C3 --> S3 + C1
    (k10, k11), S3 + C1 <--> C4
    k12, C4 --> S1 + C1
end

rn_mapk = @reaction_network begin
    (k1, k2), E + S1 <--> ES1
    k3, ES1 --> E + S2
    (k4, k5), F + S2 <--> FS2
    k6, FS2 --> F + S1
    (k7, k8), E + S2 <--> ES2
    k9, ES2 --> E + S3
    (k10, k11), F + S3 <--> FS3
    k12, FS3 --> F + S2
end

SUITE["network_construction"] = BenchmarkGroup()
SUITE["network_construction"]["simple"] = @benchmarkable @reaction_network begin
    k1, A --> B
    k2, B --> A
end
SUITE["network_construction"]["michaelis_menten"] = @benchmarkable @reaction_network begin
    (k1, k2), S + E <--> C
    k3, C --> P + E
end

SUITE["structural_analysis"] = BenchmarkGroup()
SUITE["structural_analysis"]["isconservative"] = @benchmarkable C.isconservative($rn_michaelis_menten)
SUITE["structural_analysis"]["isconsistent"] = @benchmarkable C.isconsistent($rn_michaelis_menten)
SUITE["structural_analysis"]["isconcordant"] = @benchmarkable C.isconcordant($rn_michaelis_menten)
SUITE["structural_analysis"]["minimalsiphons"] = @benchmarkable C.minimalsiphons($rn_acr)
SUITE["structural_analysis"]["elementary_flux_modes"] = @benchmarkable C.elementary_flux_modes($rn_mapk)
SUITE["structural_analysis"]["networksummary"] = @benchmarkable C.networksummary($rn_mapk)

rn_def_one = @reaction_network begin
    (k1, k2), 2A + B <--> 3A
    (k3, k4), A <--> 0
    (k5, k6), B <--> 0
end

rn_chiral = @reaction_network begin
    (k1, k2), L + 2R + P <--> 3R + Q
    (k3, k4), R + 2L + P <--> 3L + Q
    (k5, k6), P <--> 0
    (k7, k8), 0 <--> Q
end

SUITE["persistence"] = BenchmarkGroup()
SUITE["persistence"]["deficiencyone_medium"] = @benchmarkable C.deficiencyonealgorithm($rn_def_one)
SUITE["persistence"]["deficiencyone_chiral"] = @benchmarkable C.deficiencyonealgorithm($rn_chiral)
SUITE["persistence"]["ispersistent"] = @benchmarkable C.ispersistent($rn_mapk)

SUITE["concentration_robustness"] = BenchmarkGroup()
p_acr = Dict(
    [
        :k1 => 1, :k2 => 1, :k3 => 1, :k4 => 1, :k5 => 1, :k6 => 1,
        :k7 => 1, :k8 => 1, :k9 => 1, :k10 => 1, :k11 => 1, :k12 => 1,
    ]
)
SUITE["concentration_robustness"]["acr_saturation"] = @benchmarkable C.isconcentrationrobust($rn_acr)
SUITE["concentration_robustness"]["acr_with_params"] = @benchmarkable C.isconcentrationrobust($rn_acr; p = $p_acr)
