VarMapType = Union{Vector{P}, Dict, Tuple{P}} where {P <: Pair}

# Struct summarizing the dynamic information of the reaction network, including its capacity for
# multiple equilibria, concentration robustness, and persistence. 
mutable struct NetworkSummary
    steadystates::Symbol
    concentrationrobust::Symbol
    persistent::Symbol
    mixedvolume::Int
end

"""
    networksummary(rn::ReactionSystem; p::VarMapType, u0::VarMapType, mv = false)

    Summary of properties that can be inferred from the structure of a reaction network. May give different results depending on whether `p` or `u0` is supplied. The set of functions run are the following: 
    - `hasuniquesteadystates`
    - `isconcentrationrobust`
    - `ispersistent`
    - `mixedvolume`

    Mixed volume may take a very long time to run and is disabled by default. It can be enabled by setting the flag `mv = true`. Note that mixed volume requires an initial condition. 
"""
function networksummary(rn::ReactionSystem; p::VarMapType = rn.defaults, u0::VarMapType = Dict(), mv = false)
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

    NetworkSummary(eq, acr, pers, mv)
end

function Base.show(ns::NetworkSummary)
    printstyled("Number of Steady States", bold = true);
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

    println();
    printstyled("Concentration Robustness", bold = true);
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

    println();
    printstyled("Persistence", bold = true);
    println()
    if ns.persistent == :PERSISTENT
        println("This reaction network is persistent. Any species that is initially present in the reaction mixture will not die out (have its concentration reduced to zero.")
    elseif ns.persistent == :NOT_PERSISTENT
        println("This reaction network is persistent. It possess steady states that for which one or multiple species will have a concentration of zero.")
    elseif ns.persistent == :INCONCLUSIVE
        println("The algorithm currently cannot determine whether this network will have persistence.")
    end
end

"""
    hasuniquesteadystates(rn::ReactionSystem)

    Check whether a reaction network has the capacity to admit multiple steady states, for some choice of rate constants. Return codes: 
    - :NO_EQUILIBRIUM - no positive equilibrium for any choice of rate constants
    - :STRUCTURALLY_UNIQUE - only one steady state for every SCC, for every choice of rate constants
    - :STRUCTURALLY_MULTIPLE - multiple steady states in a certain SCC guaranteed for any choice of rate constants
    - :KINETICALLY_MULTIPLE - multiple steady states in a certain SCC guaranteed for certain choices of rate constants
    - :POSSIBLY_MULTIPLE - discordant and/or high deficiency, but inconclusive whether there are system parameters that lead to the existence of an SCC with multiple steady states. 
"""
function hasuniquesteadystates(rn::ReactionSystem; p::VarMapType = Dict(), u0::VarMapType = Dict())
    nps = Catalyst.get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = Catalyst.deficiency(rn)
    subs = subnetworks(rn)
    # haspositivesteadystates(rn) || return :NO_EQUILIBRIUM

    # Deficiency zero theorem 
    δ == 0 && (if isweaklyreversible(rn, subs)
        return :STRUCTURALLY_UNIQUE
    else
        return :NO_EQUILIBRIUM
    end)

    # Deficiency one networks
    Catalyst.satisfiesdeficiencyone(rn) && return :STRUCTURALLY_UNIQUE
    δ == 1 && (if deficiencyonealgorithm(rn)
        return :KINETICALLY_MULTIPLE
    else
        return :STRUCTURALLY_UNIQUE
    end)

    # Higher deficiency networks
    concordant = isconcordant(rn)
    concordant && return :STRUCTURALLY_UNIQUE
    !concordant && ispositivelydependent(rn) && return :STRUCTURALLY_MULTIPLE

    # higherdeficiencyalgorithm(rn) && return :KINETICALLY_MULTIPLE 

    # Kinetic properties
    if !isempty(p)
        length(p) != length(parameters(rn)) &&
            error("The length of the parameter map is not equal to the number of parameters in the reaction network.")
        (Catalyst.iscomplexbalanced(rn, params) ||
         Catalyst.isdetailedbalance(rn, params)) && return :KINETICALLY_UNIQUE
    end

    return :POSSIBLY_MULTIPLE
end

# Some kind of stability analysis functions?

"""
    haspositivesteadystates(rn::ReactionSystem)

    Checks whether the reaction system will have any positive steady states, i.e. steady states for which the concentration of each species is positive. 
"""
function haspositivesteadystates(rn::ReactionSystem)
    subs = subnetworks(rn)
    isweaklyreversible(rn) && return true
    !isconsistent(rn) && return false
end

"""
    haspositivesteadystates(rn::ReactionSystem)

    Checks whether the reaction system will have any periodic solutions. 
"""
function hasperiodicsolutions(rn::ReactionSystem)
    isconservative(rn) && false
end

####################################################################
# STEADY STATES IN A PARTICULAR STOICHIOMETRIC COMPATIBILITY CLASS #
####################################################################

"""
    SFR(rn::ReactionSystem; u0, p) 

    Takes in a reaction network, and returns a symbolic function that evaluates the species formation rate function, which can be used to create steady state polynomials in the desired output type (be it Symbolic, DynamicPolynomial, or QQPolyElem for abstract algebra calculations). 

    Optionally takes an initial condition (which is used to compute conservation laws) and a parameter map as arguments. These maps must be a dictionary, vector, or tuple of variable-to-value mappings, e.g. [:k1 => 1., :k2 => 2., ...]
"""
function SFR(rn::ReactionSystem; u0::VarMapType = Dict(), p::VarMapType = Dict())
    specs = species(rn);
    conslaws = conservationlaws(rn)
    sfr = if isempty(u0)
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = false, combinatoric_ratelaws = false)
    else
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = true, combinatoric_ratelaws = false)
    end

    !isa(u0, Dict) && (u0 = Dict(u0));
    !isa(p, Dict) && (p = Dict(p))

    # Substitute initial conditions. 
    if !isempty(u0)
        (length(u0) != length(specs)) &&
            error("Length of initial condition does not equal number of species.")
        u0 = symmap_to_varmap(rn, u0)
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
        p = symmap_to_varmap(rn, p)
        (length(p) != length(parameters(rn))) &&
            error("Length of parameter assignments does not equal number of parameters.")
        for i in 1:length(sfr)
            sfr[i] = Symbolics.substitute(sfr[i], p)
        end
    end

    # Generate appropriate output type. 
    argvec = vcat(species(rn), parameters(rn))
    sfr_f, sfr_f! = Symbolics.build_function(sfr, argvec...; expression = Val{false})
    sfr_f
end

"""
    modifiedSFR(rn::ReactionSystem, u0::VarMapType; p::VarMapType = Dict())

    Construct the modified SFR for the mixed volume and injectivity. Differs from the other SFR function in that certain species' rates get replaced with conservation laws, but not substituted out altogether. 
"""
function modifiedSFR(rn::ReactionSystem, u0::VarMapType; p::VarMapType = Dict())
    conslaws = conservationlaws(rn)
    d, ZZconslaws = Oscar.rref(ZZMatrix(conslaws))
    considxs = [findfirst(!=(0), conslaws[i, :]) for i in 1:d]
    conslaws = Matrix{Int64}(ZZconslaws)

    sm = speciesmap(rn);
    u0vec = zeros(length(species(rn)))
    u0 = symmap_to_varmap(rn, Dict(u0))
    for spec in keys(sm)
        i = sm[spec];
        u0vec[i] = u0[spec]
    end
    c = conslaws*u0vec

    # Get species as symbolics. 
    specs = species(rn)
    sfr = Catalyst.assemble_oderhs(rn, specs)
    conserved_eqs = conslaws*specs - c

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

    Compounds an upper bound on the number of steady states in a particular stoichiometric compatibility class. 
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
    MixedSubdivisions.mixed_volume(supp)
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
