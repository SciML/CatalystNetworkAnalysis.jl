module CatalystNetworkAnalysis

using Catalyst
using Satisfiability
using Oscar
using JuMP, HiGHS

import ModelingToolkit as MT

include("persistence.jl")
export ispersistent

end
