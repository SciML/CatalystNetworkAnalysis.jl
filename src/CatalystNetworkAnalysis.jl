module CatalystNetworkAnalysis

using Catalyst
using Satisfiability
using Oscar, Nemo
using JuMP, HiGHS
using LinearAlgebra
using Graphs
using IterTools

import ModelingToolkit as MT

const M::Float64 = 1E4
const ϵ::Float64 = 1E-1

include("persistence.jl")
export ispersistent, minimalsiphons, iscritical, isconservative, isconsistent
include("concordance.jl")
export isconcordant
include("deficiencytheory.jl")
export deficiencyonealgorithm
include("steadystates.jl")
export networksummary, modifiedSFR

end
