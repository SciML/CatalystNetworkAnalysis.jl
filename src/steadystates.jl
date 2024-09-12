
function networksummary(rn::ReactionSystem; params = rn.defaults) 
    @show rn

    # Does the network admit multiple steady states? 
    # Are any of these steady states oscillatory?
    
    # Structural Properties. 
    eq = hasuniqueequilibria(rn)
    mv = mixedvolume(rn)
    pers = ispersistent(rn)

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


    if Catalyst.satisfiesdeficiencyzero(rn)
    elseif !isempty(params)
        println("Inconclusive whether the system can admit multiple equilibria; will depend on the rate constants. The number of equilibria with this initial condition will not exceed $mv.")
    else
        println("Inconclusive whether this reaction network can admit multiple equilibria. One could try obtaining the equilibria from HomotopyContinuation.jl.") 
    end

    ### === ###
    println(); printstyled("Concentration Robustness", bold=true)
    acr = isconcentrationrobust(rn)

    if acr == :MASS_ACTION_ACR
    elseif acr == :GLOBAL_ACR
    elseif acr == :INCONCLUSIVE
    else
        println("This reaction network does not have any species that are concentration-robust.")
    end

    ### === ###
    println(); printstyled("Persistence", bold=true)
    try ispersistent(rn)
        println("This reaction network is persistent. Any species that is initially present in the reaction mixture will not die out (have its concentration reduced to zero.")
    catch error
        println("It is inconclusive whether this reaction network is persistent.")
    end
end

"""
    hasuniqueequilibria(rn::ReactionSystem)

    Check whether a reaction network has the capacity to admit multiple equilibria, for some choice of rate constants. Return codes: 
    - :NO_EQUILIBRIUM - no positive equilibrium for any choice of rate constants
    - :STRUCTURALLY_UNIQUE - only one steady state for every SCC, for every choice of rate constants
    - :KINETICALLY_UNIQUE - only one steady state for every SCC, for this particular set of rate constants
    - :DEFINITELY_MULTIPLE - multiple steady states in a certain SCC guaranteed for any choice of rate constants
    - :KINETICALLY_MULTIPLE - multiple steady states in a certain SCC guaranteed for certain choices of rate constants
    - :POSSIBLY_MULTIPLE - discordant and/or high deficiency, but inconclusive whether there are system parameters that lead to the existence of an SCC with multiple steady states. 
"""

function hasuniqueequilibria(rn::ReactionSystem, params) 
    nps = get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    haspositiveequilibria(rn) || error("This reaction network does not have the ability to admit positive equilibria for any choice of rate constants.")
    haspositiveequilibria(rn) || return :NO_EQUILIBRIUM

    # Deficiency zero theorem 
    δ == 0 && (isweaklyreversible(rn) ? return :STRUCTURALLY_UNIQUE : return :NO_EQUILIBRIUM)

    # Deficiency one networks
    Catalyst.satisfiesdeficiencyone(rn) && return :STRUCTURALLY_UNIQUE 
    δ == 1 && (deficiencyonealgorithm(rn) ? return :KINETICALLY_MULTIPLE : return :STRUCTURALLY_UNIQUE)

    # Higher deficiency networks
    concordant = isconcordant(rn)
    concordant && return :STRUCTURALLY_UNIQUE 
    !concordant && ispositivelydependent(rn) && return :DEFINITELY_MULTIPLE 

    # Kinetic properties
    (Catalyst.iscomplexbalanced(rn, params) || Catalyst.isdetailedbalance(rn, params)) && return :KINETICALLY_UNIQUE  
    higherdeficiencyalgorithm(rn) && return :KINETICALLY_MULTIPLE 
    
    return :POSSIBLY_MULTIPLE
    error("The network is discordant and high deficiency, but this function currently cannot conclude whether the network has the potential to have multiple equilibria.")
end

"""
    isconcentrationrobust(rn::ReactionSystem)

    Check whether a reaction network has any concentration-robust species. Return codes: 
    - :MASS_ACTION_ACR - this species is concentraiton-robust for the given set of rate constants 
    - :UNCONDITIONAL_ACR - this species is absolutely concentraiton-robust for every choice of rate constants
    - :INCONCLUSIVE - the algorithm currently cannot decide whether this network has ACR. One could try calling this function with rate constants provided. 
"""

function isconcentrationrobust(rn::ReactionSystem) 
    
end

function hasperiodicsolutions(rn::ReactionSystem) 
    isconservative(rn) && false

    error("Inconclusive.")
end

# Some kind of stability analysis functions?
function haspositiveequilibria(rn::ReactionSystem) 
    
end
 
# Upper bound on the number of steady states. TODO: get this to work with polynomials in Symbolics
function mixedvolume(rn::ReactionSystem) 
    conslaws = conservationlaws(rn); (d, s) = size(conslaws)
    idxs = Set([findfirst(!=(0), conslaws[i, :]) for i in 1:d])

    specs = species(rn)
    Wx_c = conslaws*specs
end
