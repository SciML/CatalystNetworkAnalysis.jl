module CatalystNetworkAnalysis

using Catalyst
using Satisfiability
using Oscar, Nemo
using JuMP, HiGHS
using LinearAlgebra
using Graphs
using IterTools

import ModelingToolkit as MT

include("persistence.jl")
export ispersistent, minimalsiphons, iscritical, isconservative, isconsistent
include("concordance.jl")
export isconcordant
include("deficiencyonealgorithm.jl")
export deficiencyonealgorithm

end
