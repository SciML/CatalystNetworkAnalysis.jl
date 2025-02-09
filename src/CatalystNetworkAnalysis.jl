module CatalystNetworkAnalysis

using PrecompileTools: @setup_workload, @compile_workload

using Catalyst
using Satisfiability # For siphon detection

# Algebraic functionality
using Oscar

using JuMP, HiGHS # For concordance and deficiency algorithms
const M::Float64 = 1E6
const ϵ::Float64 = 1E-6 # Constants for the linear programming solvers 

using LinearAlgebra
using Graphs
using IterTools, Combinatorics
using SparseArrays

using MixedSubdivisions, DynamicPolynomials # For polytope analysis
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
include("steady_state_parameterizations.jl")
export symbolic_steady_states

end
