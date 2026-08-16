VarMapType = Union{Vector{P}, Dict, Tuple{P}} where {P <: Pair}

# Catalyst ≥ 16 (ModelingToolkitBase) replaced the system `defaults` field with
# `initial_conditions`; older Catalyst exposes the defaults dict via `rn.defaults`.
@static if isdefined(Catalyst, :initial_conditions)
    systemdefaults(rn::ReactionSystem) = Dict(Catalyst.initial_conditions(rn))
else
    systemdefaults(rn::ReactionSystem) = rn.defaults
end

"""
    NetworkSummary

Summary of structural properties computed by [`networksummary`](@ref).

# Fields

- `steadystates`: status returned by `hasuniquesteadystates`.
- `concentrationrobust`: status returned by `isconcentrationrobust`.
- `persistent`: status returned by `ispersistent`.
- `mixedvolume`: upper bound from `mixedvolume`, or `-1` when that calculation
  was not requested.
"""
mutable struct NetworkSummary
    steadystates::Symbol
    concentrationrobust::Symbol
    persistent::Symbol
    mixedvolume::Int
end

"""
    networksummary(rn::ReactionSystem; p::VarMapType, u0::VarMapType, mv = false)

Compute a structural summary of a reaction system.

# Arguments

- `rn`: mass-action reaction system to analyze.

# Keyword Arguments

- `p`: parameter map passed to the steady-state and concentration-robustness
  analyses. The system defaults are used by default.
- `u0`: initial-condition map used by the mixed-volume analysis.
- `mv`: whether to compute the mixed-volume bound. This can be expensive and
  defaults to `false`; it is ignored when `u0` is empty.

# Returns

A [`NetworkSummary`](@ref) containing the steady-state, concentration-
robustness, persistence, and mixed-volume statuses.

# Throws

An error if `rn` is not an integer-coefficient mass-action reaction system.
"""
function networksummary(rn::ReactionSystem; p::VarMapType = systemdefaults(rn), u0::VarMapType = Dict(), mv = false)
    all(r -> ismassaction(r, rn), reactions(rn)) ||
        error("The network summary analysis currently only works for mass-action networks with integer coefficients.")

    # Structural Properties.
    eq = hasuniquesteadystates(rn; p = p)
    acr = isconcentrationrobust(rn; p = p)
    mv = if (isempty(u0) || !mv)
        -1
    else
        mixedvolume(rn, u0)
    end
    pers = ispersistent(rn)

    return NetworkSummary(eq, acr, pers, mv)
end

function Base.show(ns::NetworkSummary)
    printstyled("Number of Steady States", bold = true)
    println()
    if ns.steadystates == :STRUCTURALLY_UNIQUE
        println("This reaction network will have a unique steady-state for every stoichiometric compatibility class, for every choice of rate constants. If the network is deficiency zero, this steady-state will additionally be asymptotically stable.")
    elseif ns.steadystates == :KINETICALLY_UNIQUE
        println("The number of steady states will depend on the rate constants. For the choice given,")
    elseif ns.steadystates == :POSSIBLY_MULTIPLE
        println("Inconclusive whether the system can admit multiple steady states; will depend on the rate constants. One could try obtaining the steady states from HomotopyContinuation. Note that the number of steady states depends in general on the choice of initial condition.")
        if ns.mv >= 0
            println("The number of steady states with this initial condition will not exceed $(ns.mv).")
        end
    elseif ns.steadystates == :DEFINITELY_MULTIPLE
        println("This network is guaranteed to have an initial condition for which there are multiple steady states, for any choice of rate constants. Try running a stability analysis.")
    elseif ns.steadystates == :KINETICALLY_MULTIPLE
        println("This network is guaranteed to have an initial condition for which there are multiple steady states, for a certain set of rate constants.") # Try running ? for an example set of rate constants.

    elseif ns.steadystates == :NO_EQUILIBRIUM
        println("This reaction network will not have positive steady states, for any choice of rate constants.")
    else
        error("Unrecognized status message for multiple steady states.")
    end

    println()
    printstyled("Concentration Robustness", bold = true)
    println()

    if ns.concentrationrobust == :MASS_ACTION_ACR
        println("This reaction network has absolute concentration robustness in at least one species for this set of rate parameters. The concentration of this species will be constant ns.concentrationrobustoss all steady states for the system. To see the indices of species that are concentration robust, please query nps.robustspecies.")
    elseif ns.concentrationrobust == :GLOBAL_ACR
        println("This reaction network has absolute concentration robustness in at least one species for any set of rate parameters. The concentration of this species will be constant ns.concentrationrobustoss all steady states for the system. To see the indices of species that are concentration robust, please query nps.robustspecies.")
    elseif ns.concentrationrobust == :INCONCLUSIVE
        println("The algorithm currently cannot determine whether this network will have concentration robustness in any species.")
    elseif ns.concentrationrobust == :NO_ACR
        println("This reaction network does not possess absolute concentration robustness in any species, for any set of rate constants.")
    end

    println()
    printstyled("Persistence", bold = true)
    println()
    return if ns.persistent == :PERSISTENT
        println("This reaction network is persistent. Any species that is initially present in the reaction mixture will not die out (have its concentration reduced to zero.")
    elseif ns.persistent == :NOT_PERSISTENT
        println("This reaction network is persistent. It possess steady states that for which one or multiple species will have a concentration of zero.")
    elseif ns.persistent == :INCONCLUSIVE
        println("The algorithm currently cannot determine whether this network will have persistence.")
    end
end

"""
    hasuniquesteadystates(rn::ReactionSystem)

Classify the possible number of positive steady states of a reaction system.

# Arguments

- `rn`: reaction system to analyze.

# Keyword Arguments

- `p`: optional parameter map used for kinetic balance checks.
- `u0`: optional initial-condition map reserved for compatibility-class
  analyses.

# Returns

One of `:NO_EQUILIBRIUM`, `:STRUCTURALLY_UNIQUE`,
`:STRUCTURALLY_MULTIPLE`, `:KINETICALLY_UNIQUE`, `:KINETICALLY_MULTIPLE`, or
`:POSSIBLY_MULTIPLE`.
"""
function hasuniquesteadystates(rn::ReactionSystem; p::VarMapType = Dict(), u0::VarMapType = Dict())
    nps = Catalyst.get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = Catalyst.deficiency(rn)
    subs = subnetworks(rn)
    # haspositivesteadystates(rn) || return :NO_EQUILIBRIUM

    # Deficiency zero theorem
    δ == 0 && (
        if isweaklyreversible(rn, subs)
            return :STRUCTURALLY_UNIQUE
        else
            return :NO_EQUILIBRIUM
        end
    )

    # Deficiency one networks
    Catalyst.satisfiesdeficiencyone(rn) && return :STRUCTURALLY_UNIQUE
    δ == 1 && (
        if deficiencyonealgorithm(rn)
            return :KINETICALLY_MULTIPLE
        else
            return :STRUCTURALLY_UNIQUE
        end
    )

    # Higher deficiency networks
    concordant = isconcordant(rn)
    concordant && return :STRUCTURALLY_UNIQUE
    !concordant && ispositivelydependent(rn) && return :STRUCTURALLY_MULTIPLE

    # higherdeficiencyalgorithm(rn) && return :KINETICALLY_MULTIPLE

    # Kinetic properties
    if !isempty(p)
        length(p) != length(parameters(rn)) &&
            error("The length of the parameter map is not equal to the number of parameters in the reaction network.")
        (
            Catalyst.iscomplexbalanced(rn, params) ||
                Catalyst.isdetailedbalance(rn, params)
        ) && return :KINETICALLY_UNIQUE
    end

    return :POSSIBLY_MULTIPLE
end

# Some kind of stability analysis functions?

"""
    haspositivesteadystates(rn::ReactionSystem)

Check whether a reaction system admits a positive steady state.

# Arguments

- `rn`: reaction system to analyze.

# Returns

A Boolean result when the structural checks decide the question. The current
implementation is conservative for cases that require a more detailed
analysis.
"""
function haspositivesteadystates(rn::ReactionSystem)
    subs = subnetworks(rn)
    isweaklyreversible(rn) && return true
    return !isconsistent(rn) && return false
end

"""
    hasperiodicsolutions(rn::ReactionSystem)

Check whether a reaction system admits periodic solutions.

# Arguments

- `rn`: reaction system to analyze.

# Returns

`false` until a periodic-solution criterion is implemented. This developer
hook is retained so callers can use a stable predicate while the analysis is
expanded.
"""
function hasperiodicsolutions(rn::ReactionSystem)
    return isconservative(rn) && false
end

####################################################################
# STEADY STATES IN A PARTICULAR STOICHIOMETRIC COMPATIBILITY CLASS #
####################################################################

"""
    SFR(rn::ReactionSystem; u0 = Dict(), p = Dict())

Construct a callable species-formation-rate function.

The returned function accepts the species followed by the parameters of `rn`.
It can therefore be evaluated in the symbolic or polynomial type required by
downstream steady-state analyses.

# Arguments

- `rn`: reaction system whose formation rates are assembled.

# Keyword Arguments

- `u0`: optional initial-condition map used to substitute conservation-law
  constants.
- `p`: optional parameter map used to substitute rate parameters.

Both maps may be a `Dict`, a vector of `Pair`, or a tuple of `Pair`.

# Returns

A callable function of the species and parameters in `rn`.
"""
function SFR(rn::ReactionSystem; u0::VarMapType = Dict(), p::VarMapType = Dict())
    specs = species(rn)
    conslaws = conservationlaws(rn)
    sfr = if isempty(u0)
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = false, combinatoric_ratelaws = false)
    else
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = true, combinatoric_ratelaws = false)
    end

    !isa(u0, Dict) && (u0 = Dict(u0))
    !isa(p, Dict) && (p = Dict(p))

    # Substitute initial conditions.
    if !isempty(u0)
        (length(u0) != length(specs)) &&
            error("Length of initial condition does not equal number of species.")
        u0 = Catalyst.symmap_to_varmap(rn, u0)
        cons_constants = Catalyst.conservationlaw_constants(rn)
        Γ_vals = Vector{Float64}()
        for conseq in cons_constants
            push!(Γ_vals, Symbolics.substitute(conseq.rhs, u0))
        end
        cons_map = Dict(cons.lhs => Γ_val for (cons, Γ_val) in zip(cons_constants, Γ_vals))
        for i in 1:length(sfr)
            sfr[i] = Symbolics.substitute(sfr[i], cons_map)
        end
    end

    # Substitute parameters.
    if !isempty(p)
        p = Catalyst.symmap_to_varmap(rn, p)
        (length(p) != length(parameters(rn))) &&
            error("Length of parameter assignments does not equal number of parameters.")
        for i in 1:length(sfr)
            sfr[i] = Symbolics.substitute(sfr[i], p)
        end
    end

    # Generate appropriate output type.
    argvec = vcat(species(rn), parameters(rn))
    sfr_f, sfr_f! = Symbolics.build_function(sfr, argvec...; expression = Val{false})
    return sfr_f
end

"""
    modifiedSFR(rn::ReactionSystem, u0::VarMapType; p::VarMapType = Dict())

Construct the modified species-formation-rate system used by mixed-volume and
injectivity analyses.

Unlike [`SFR`](@ref), this representation replaces selected species rates by
conservation-law equations without eliminating those species from the system.

# Arguments

- `rn`: reaction system to analyze.
- `u0`: initial-condition map used to evaluate the conserved quantities.

# Keyword Arguments

- `p`: parameter map reserved for the parameterized formation-rate workflow.

# Returns

The modified symbolic formation-rate vector with the conservation equations
appended for the selected species.
"""
function modifiedSFR(rn::ReactionSystem, u0::VarMapType; p::VarMapType = Dict())
    conslaws = conservationlaws(rn)
    d, ZZconslaws = Oscar.rref(ZZMatrix(conslaws))
    considxs = [findfirst(!=(0), conslaws[i, :]) for i in 1:d]
    conslaws = Matrix{Int64}(ZZconslaws)

    sm = speciesmap(rn)
    u0vec = zeros(length(species(rn)))
    u0 = Catalyst.symmap_to_varmap(rn, Dict(u0))
    for spec in keys(sm)
        i = sm[spec]
        u0vec[i] = u0[spec]
    end
    c = conslaws * u0vec

    # Get species as symbolics.
    specs = species(rn)
    sfr = Catalyst.assemble_oderhs(rn, specs)
    conserved_eqs = conslaws * specs - c

    for (i, rx) in enumerate(considxs)
        sfr[rx] = conserved_eqs[i]
    end

    argvec = vcat(species(rn), parameters(rn))
    sfr_f, sfr_f! = Symbolics.build_function(sfr, argvec...; expression = Val{false})
    sfr_f
    return sfr_f
end

"""
    mixedvolume(rn::ReactionSystem, u0::VarMapType)

Compute an upper bound on the number of steady states in a compatibility class.

# Arguments

- `rn`: reaction system to analyze.
- `u0`: initial-condition map defining the stoichiometric compatibility class.

# Returns

The mixed volume of the modified steady-state polynomial support.

# Throws

An error if `u0` does not contain one value for every species in `rn`.
"""
function mixedvolume(rn::ReactionSystem, u0::VarMapType)
    (length(u0) != length(species(rn))) &&
        error("The length of the initial condition must equal the number of species in the reaction network.")
    sfr_f = modifiedSFR(rn, u0)
    sfr = eval(sfr_f)

    @polyvar s[1:length(species(rn))]
    @polyvar k[1:length(parameters(rn))]

    polysfr = sfr(vcat(s, k)...)
    supp = MixedSubdivisions.support(polysfr, s)
    return MixedSubdivisions.mixed_volume(supp)
end

# """
#     isinjective(rn, u0)
# """
# function isinjective(rn::ReactionSystem, u0::Vector)
#     # check ispermanent(rn)
#     sfr = modifiedSFR(rn, u0)
#     J = Symbolics.jacobian(sfr, species(rn))
#     detJ = det(J)
#
#     # Check positivity.
# end
#
# # Steady States in an SCC.
