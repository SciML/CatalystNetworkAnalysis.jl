
function isconcentrationrobust(rn::ReactionNetwork) 
    nps = get_networkproperties(rn)
    # Compute steady state ideal
    R, specvars = polynomial_ring(QQ, map(s -> Symbolics.tosymbol(s, escape=false), species(rn)))
    sfr = symbolicSFR(rn, specvars)
    ss_ideal = 0  

    # If Grobner basis has an element of the form (x_i - α), we have ACR in i
    I = ideal(R, sfr)
    le = linearelements(I)
    if !isempty(le)
        nps.isrobust = true
        push!(nps.robustspecies, le...)
        return true
    end

    # Compute saturation of ideal
    J = prod(specvars)
    satI = saturation(I, J)
    le = linearelements(satI)
    if !isempty(le)
        nps.isrobust = true
        push!(nps.robustspecies, le...)
        return true
    end
    # Iterate over elimination ideals I ∩ Q[x_i] and check for unique roots
    
    # Compute the steady-state parameterization, check for constant components. 
    
    # Compute decompositions of the positive-restriction ideal
    
end

function linearelements(I::Ideal) 
    G = grobner_basis(I)
    linearidxs = Vector{Int64}()
    for g in elements(G)
        deg = degrees(g)
        islinear = (count(==(0), deg) == length(deg) - 1 && count(==(1), deg) == 1)
        if islinear
            push!(linearidxs, findfirst(==(1), deg))
        end
    end
    return linearidxs
end

function robustspecies(rn::ReactionNetwork) 

end    

function positiverestrictionideal(I::Ideal{T}) where T <: MPolyRingElem

end
