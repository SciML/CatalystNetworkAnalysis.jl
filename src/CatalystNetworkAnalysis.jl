module CatalystNetworkAnalysis

using Catalyst
using Satisfiability # For siphon detection
using Oscar, Nemo, Hecke # Algebraic functionality
using JuMP, HiGHS # For concordance and deficiency algorithms
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
export networksummary, SFR
include("concentrationrobustness.jl")
export isconcentrationrobust

end
