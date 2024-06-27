module CatalystNetworkAnalysis

using Catalyst
using Satisfiability
using Oscar
using JuMP, HiGHS
using LinearAlgebra

import ModelingToolkit as MT

include("persistence.jl")
export ispersistent, minimalsiphons
include("concordance.jl")
export isconcordant
include("deficiencyonealgorithm.jl")
export deficiencyonealgorithm


end
