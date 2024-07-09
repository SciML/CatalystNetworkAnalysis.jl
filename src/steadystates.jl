"""
    hasuniqueequilibria(rn::ReactionSystem)

    Check whether a reaction network has the capacity to admit multiple equilibria, for some choice of rate constants. 
"""

function hasuniqueequilibria(rn::ReactionSystem) 
    nps = get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    concordant = isconcordant(rn)

    δ == 0 && return true 
    Catalyst.satisfiesdeficiencyone(rn) && return true 
    concordant && return true
    !concordant && ispositivelydependent(rn) && return false
    δ == 1 && return deficiencyonealgorithm(rn)
    
    error("The network is discordant and high deficiency, but this function currently cannot conclude whether the network has the potential to have multiple equilibria.")
end


function hasperiodicsolutions(rn::ReactionSystem) 
    
end

# Some kind of stability analysis functions?
# 
