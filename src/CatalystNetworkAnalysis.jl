module CatalystNetworkAnalysis

using Catalyst
using PrecompileTools: @compile_workload, @setup_workload
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

import SymbolicUtils

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
export WRDZ_translation

@setup_workload begin
    @compile_workload begin
        rn = @reaction_network begin
            k₁, A --> B
            k₂, B --> A
        end
        isconservative(rn)
        isconsistent(rn)
    end
end

end
