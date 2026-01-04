module CatalystNetworkAnalysis

using PrecompileTools: @setup_workload, @compile_workload

using Catalyst
using Satisfiability # For siphon detection

# Algebraic functionality
using Oscar

# Linear programming (for concordance + deficiency)
using JuMP, HiGHS
const M::Float64 = 1.0e4
const ϵ::Float64 = 1.0e-4

using LinearAlgebra
using Graphs
using IterTools, Combinatorics
using SparseArrays

# Polytope analysis (EFMs)
using MixedSubdivisions, DynamicPolynomials
using Polyhedra
import CDDLib

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
