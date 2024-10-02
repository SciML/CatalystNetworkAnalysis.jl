
function isconcentrationrobust(rn::ReactionNetwork) 
    # Compute steady state ideal
    R, specvars = polynomial_ring(QQ, map(s -> Symbolics.tosymbol(s, escape=false), species(rn)))
    sfr = symbolicSFR(rn, specvars)
    ss_ideal = 0  
    # If Grobner basis has an element of the form (x_i - α), we have ACR in i
    
    # Compute saturation of ideal
    # Iterate over elimination ideals I ∩ Q[x_i] and check for unique roots
    
    # Compute the steady-state parameterization, check for constant components. 
    
    # Compute decompositions of the positive-restriction ideal
    
end

function robustspecies(rn::ReactionNetwork) 

end    

function positiverestrictionideal(I::Ideal{T}) where T <: MPolyRingElem

end
