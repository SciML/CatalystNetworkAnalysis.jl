"""
    hasuniqueequilibria(rn::ReactionSystem)

    Check whether a reaction network has the capacity to admit multiple equilibria, for some choice of rate constants. 
"""

function hasuniqueequilibria(rn::ReactionSystem) 
    nps = get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    haspositiveequilibria(rn) || error("This reaction network does not have the ability to admit positive equilibria for any choice of rate constants.")

    # Deficiency one theorem 
    δ == 0 && return true 

    # Deficiency one networks
    Catalyst.satisfiesdeficiencyone(rn) && return true 
    δ == 1 && return deficiencyonealgorithm(rn)

    # Higher deficiency networks
    concordant = isconcordant(rn)
    concordant && return true
    !concordant && ispositivelydependent(rn) && return false
    
    error("The network is discordant and high deficiency, but this function currently cannot conclude whether the network has the potential to have multiple equilibria.")
end


function hasperiodicsolutions(rn::ReactionSystem) 
    
end

# Some kind of stability analysis functions?
function haspositiveequilibria(rn::ReactionSystem) 
    
end
 
