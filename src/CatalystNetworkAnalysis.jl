module CatalystNetworkAnalysis

using Catalyst
using Satisfiability # For siphon detection

# Algebraic functionality
using Oscar, Nemo
import Hecke: n_positive_roots

using JuMP, HiGHS # For concordance and deficiency algorithms
using MixedSubdivisions, DynamicPolynomials # For polytope analysis
using LinearAlgebra
using Graphs

using SparseArrays, StaticArrays
using IterTools, Combinatorics

using Polyhedra
import CDDLib

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

include("utils.jl")
include("lp_utils.jl")
include("cycles.jl")
export elementary_flux_modes

include("translated.jl") 
export WRDZ_translation, symbolic_steady_states

end
