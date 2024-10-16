# Struct summarizing the dynamic information of the reaction network, including its capacity for
# multiple equilibria, concentration robustness, and persistence. 

mutable struct NetworkSummary
    Equilibria::Enum
    ConcentrationRobust::Enum
    Persistent::Enum
end

function networksummary(rn::ReactionSystem; params = rn.defaults) 
    # Does the network admit multiple steady states? 
    # Are any of these steady states oscillatory?
    
    # Structural Properties. 
    eq = hasuniquesteadystates(rn)
    acr = isconcentrationrobust(rn)
    # mv = mixedvolume(rn)
    pers = ispersistent(rn)

    NetworkSummary(eq, acr, pers)
end

function Base.show(ns::NetworkSummary) 
    printstyled("Number of Steady States", bold=true)
    if eq == :STRUCTURALLY_UNIQUE
        println("This reaction network will have a unique steady-state for every stoichiometric compatibility class, for every choice of rate constants. If the network is deficiency zero, this steady-state will additionally be asymptotically stable.")
    elseif eq == :KINETICALLY_UNIQUE
        println("The number of steady states will depend on the rate constants. For the choice given,")
    elseif eq == :POSSIBLY_MULTIPLE
        println("Inconclusive whether the system can admit multiple steady states; will depend on the rate constants. The number of equilibria with this initial condition will not exceed $mv. One could try obtaining the steady states from HomotopyContinuation. Note that the number of steady states depends in general on the choice of initial condition.") 
    elseif eq == :DEFINITELY_MULTIPLE
        println("This network is guaranteed to have an initial condition for which there are multiple steady states, for any choice of rate constants. Try running a stability analysis.")
    elseif eq == :KINETICALLY_MULTIPLE
        println("This network is guaranteed to have an initial condition for which there are multiple steady states, for a certain set of rate constants.") # Try running ? for an example set of rate constants. 

    elseif eq == :NO_EQUILIBRIUM
        println("This reaction network will not have positive steady states, for any choice of rate constants.")
    else
        error("Unrecognized status message for multiple equilibria.")
    end

    println(); printstyled("Concentration Robustness", bold=true)
    acr = isconcentrationrobust(rn)

    if acr == :MASS_ACTION_ACR
    elseif acr == :GLOBAL_ACR
    elseif acr == :INCONCLUSIVE
    else
        println("This reaction network does not have any species that are concentration-robust.")
    end

    println(); printstyled("Persistence", bold=true)
    try ispersistent(rn)
        println("This reaction network is persistent. Any species that is initially present in the reaction mixture will not die out (have its concentration reduced to zero.")
    catch error
        println("It is inconclusive whether this reaction network is persistent.")
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

function hasuniquesteadystates(rn::ReactionSystem, params) 
    nps = get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    # haspositivesteadystates(rn) || error("This reaction network does not have the ability to admit positive equilibria for any choice of rate constants.")
    haspositivesteadystates(rn) || return :NO_EQUILIBRIUM

    # Deficiency zero theorem 
    δ == 0 && (if isweaklyreversible(rn) 
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

    # Kinetic properties
    # (Catalyst.iscomplexbalanced(rn, params) || Catalyst.isdetailedbalance(rn, params)) && return :KINETICALLY_UNIQUE  
    # higherdeficiencyalgorithm(rn) && return :KINETICALLY_MULTIPLE 
    
    return :POSSIBLY_MULTIPLE
end

"""
    isconcentrationrobust(rn::ReactionSystem)

    Check whether a reaction network has any concentration-robust species. Return codes: 
    - :MASS_ACTION_ACR - this species is concentration-robust for the given set of rate constants 
    - :UNCONDITIONAL_ACR - this species is absolutely concentration-robust for every choice of rate constants
    - :INCONCLUSIVE - the algorithm currently cannot decide whether this network has ACR. One could try calling this function with rate constants provided. 
"""


# Some kind of stability analysis functions?
function haspositivesteadystates(rn::ReactionSystem) 
    isweaklyreversible(rn) && return true
end
 
# Check whether a reaction network has periodic solutions. 
function hasperiodicsolutions(rn::ReactionSystem) 
    isconservative(rn) && false

    error("Inconclusive.")
end

####################################################################
# STEADY STATES IN A PARTICULAR STOICHIOMETRIC COMPATIBILITY CLASS #
####################################################################


# Need parameters for the symbolic SFR. 
function symbolicSFR(rn::ReactionSystem, vars::Vector{T}) where T <: MPolyRingElem
    rxs = reactions(rn)
    !all(rx -> ismassaction(rx, rn), rxs) && 
        error("Symbolic species-formation-rate in terms of polynomial ring elements is currently only supported for mass action systems.")

    spectoidx = Dict(x => i for (i, x) in enumerate(species(rn)))
    sfr = vars - vars

    for rx in rxs
        rl = rx.rate
        for (spec, stoich) in rx.substoich
        end
        for (spec, stoich) in rx.netstoich
            i = spectoidx[spec]

        end
    end
end

"""
    SFR(rn::ReactionSystem) 
    
    Takes in a reaction network, and returns a symbolic function of the species formation rate function,
    which can be used to create steady state polynomials in the desired output type (be it Symbolic, DynamicPolynomial, or QQPolyElem for abstract algebra calculations). 
"""
function SFR(rn::ReactionSystem, u0::Dict = Dict(), p::Dict = Dict(); output = :SYM) 
    specs = species(rn); conslaws = conservationlaws(rn) 
    sfr = if isempty(u0)
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = false)
    else
        Catalyst.assemble_oderhs(rn, specs, remove_conserved = true)
    end
    
    # Substitute initial conditions. 
    if !isempty(u0) 
        (length(u0) != length(specs)) && error("Length of initial condition does not equal number of species.")
        u0 = symmap_to_varmap(rn, u0)
        cons_constants = Catalyst.conservationlaw_constants(rn)
        Γ_vals = Vector{Float64}()
        for conseq in cons_constants
            push!(Γ_vals, substitute(conseq.rhs, u0))
        end
        cons_map = Dict(cons.lhs => Γ_val for (cons, Γ_val) in zip(cons_constants, Γ_vals))
        for i in 1:length(sfr)
            sfr[i] = substitute(sfr[i], cons_map)
        end
    end

    # Substitute parameters.
    if !isempty(p)
        p = symmap_to_varmap(rn, p)
        (length(p) != length(parameters(rn))) && error("Length of parameter assignments does not equal number of parameters.")
        for i in 1:length(sfr)
            sfr[i] = substitute(sfr[i], p)
        end
    end

    # Generate appropriate output type. 
    argvec = vcat(species(rn), parameters(rn))
    sfr_f, sfr_f! = Symbolics.build_function(sfr, argvec...)
    sfr_f
end

function SFR(rn::ReactionSystem, u0::Vector, p::Vector) 
    u0map = Dict([spec => u0_i for (spec, u0_i) in zip(species(rn), u0)])
    pmap = Dict([param => p_i for (param, p_i) in zip(parameters(rn), p)])
    SFR_expr(rn, u0map, pmap)
end

"""
    Macro that evaluates the SFR expression, using variables of the desired type (DP, Oscar, Symbolics, etc.)
"""

function modifiedSFR(rn::ReactionSystem, u0::Vector{Float64}) 
    conslaws = conservationlaws(rn) 
    d, ZZconslaws = Oscar.rref(ZZMatrix(conslaws))
    considxs = [findfirst(!=(0), conslaws[i, :]) for i in 1:d]
    conslaws = Matrix{Int64}(ZZconslaws)
    c = conslaws*u0

    # Get species as symbolics. 
    specs = species(rn)
    sfr = Catalyst.assemble_oderhs(rn, specs)
    conserved_eqs = conslaws*specs - c
    
    for (i, rx) in enumerate(considxs)
        sfr[rx] = conserved_eqs[i]
    end
    return sfr
end

# Upper bound on the number of steady states in a particular stoichiometric compatibility class. 
# TODO: test if this works with non-mass action kinetics. Should error if the rate constants contain function terms. 
function mixedvolume(rn::ReactionSystem, u0) 
    sfr = modifiedSFR(rn, u0)
    pvar2sym, sym2term = SymbolicUtils.get_pvar2sym(), SymbolicUtils.get_sym2term()    
    polysfr = map(eq -> PolyForm(eq, pvar2sym, sym2term).p, sfr)

    # Compute mixed volume. Get the species terms. 
    specvarnames = collect(keys(sym2term))
    vars = [pvar2sym(specvar) for specvar in specvarnames]
    supp = support(polysfr, vars)
    mixed_volume(supp)
end

# TODO
function isinjective(rn::ReactionSystem, u0::Vector) 
    # check ispermanent(rn)
    sfr = modifiedSFR(rn, u0)
    J = Symbolics.jacobian(sfr, species(rn))
    detJ = det(J)

    # Check positivity. 
end

# Steady States in an SCC. 
