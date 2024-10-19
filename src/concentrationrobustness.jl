"""
    isconcentrationrobust(rn::ReactionSystem)

    Requires the parameter values be Rational or Integer. Check whether a reaction network has any concentration-robust species. Return codes: 
    - :MASS_ACTION_ACR - this species is concentraiton-robust for the given set of rate constants 
    - :GLOBAL_ACR - this species is absolutely concentraiton-robust for every choice of rate constants
    - :INCONCLUSIVE - the algorithm currently cannot decide whether this network has ACR. One could try calling this function with rate constants provided. 
    - :NO_ACR - the reaction network does not have ACR. 
"""

function isconcentrationrobust(rn::ReactionSystem; p::Dict{Any, Rational} = Dict()) 
    nps = get_networkproperties(rn)

    # Need RATIONAL values for parameters. 
    for param in keys(p)
        p[param] = rationalize(p[param])
    end

    # Compute steady state ideal. 
    sfr_f = eval(C.SFR(rn, Dict(), p))
    R, polyvars = polynomial_ring(QQ, vcat(
                     map(s -> Symbolics.tosymbol(s, escape=false), species(rn)), 
                     map(p -> Symbolics.tosymbol(p), parameters(rn))
                 )) 
    sfr = sfr_f(polyvars...)
    specs = polyvars[1:numspecies(rn)]

    # If Grobner basis has an element of the form (x_i - α), we have ACR in i
    I = ideal(R, sfr)
    le = linearelements(I, length(specs))
    if !isempty(le)
        nps.isrobust = true
        push!(nps.robustspecies, le...)
        return isempty(p) ? :GLOBAL_ACR : :MASS_ACTION_ACR
    end

    # Compute saturation of ideal
    J = prod(specs)
    satI = saturation(I, J)
    le = linearelements(satI, length(specs))
    if !isempty(le)
        nps.isrobust = true
        push!(nps.robustspecies, le...)
        return isempty(p) ? :GLOBAL_ACR : :MASS_ACTION_ACR
    end

    # Iterate over elimination ideals I ∩ Q[x_i] and check for unique positive root. TODO: Need to write a check for unique positive roots. 

    # Create a dummy univariate ring
    r_ξ, ξ = polynomial_ring(QQ, :ξ)

    !isempty(p) && for i in 1:numspecies(rn)
        IQ = eliminate(I, vcat(specs[1:i-1], specs[i+1:end]))
        for g in gens(IQ)
            coeffs = [exponent_vector(g, j)[i] for j in 1:length(Oscar.terms(g))]
            poly = r_ξ(coeffs)
            iszero(poly) && continue 
            if Hecke.n_positive_roots(poly) == 1
                push!(nps.robustspecies, i)
                return :MASS_ACTION_ACR
            end
        end
    end
    
    # Compute decompositions of the positive-restriction ideal

    # Numerically checking for ACR
    
    # Compute the steady-state parameterization, check for constant components. 
    
    # Algorithm 6.1 for finding all pairs of ACR candidates
    
    return :INCONCLUSIVE
end

"""
    linearelements(I, numspecies)

    Look for terms of the form x - α_i in the basis for the positive steady-state ideal. Their presence
    implies the existence of ACR in that species. 
"""
function linearelements(I::Ideal, numspecies::Int) 
    G = groebner_basis(I)
    linearidxs = Vector{Int64}()
    for g in elements(G)
        deg = degrees(g)[1:numspecies]
        islinear = (count(==(0), deg) == length(deg) - 1 && count(==(1), deg) == 1)
        if islinear
            push!(linearidxs, findfirst(==(1), deg))
        end
    end
    return linearidxs
end


# Other methods. 

function positiverestrictionideal(I::Ideal{T}) where T <: MPolyRingElem

end

function denseSampling() 
    
end
