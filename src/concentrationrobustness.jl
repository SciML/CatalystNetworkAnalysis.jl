"""
    isconcentrationrobust(rn::ReactionSystem)

    Requires the parameter values be Rational or Integer. Check whether a reaction network has any concentration-robust species. Return codes: 
    - :MASS_ACTION_ACR - this species is concentration-robust for the given set of rate constants 
    - :GLOBAL_ACR - this species is absolutely concentration-robust for every choice of rate constants
    - :INCONCLUSIVE - the algorithm currently cannot decide whether this network has ACR. One could try calling this function with rate constants provided. 
    - :NO_ACR - the reaction network does not have ACR. 
    - :INEXACTPARAMS - the algorithm cannot conclude concentration-robustness due to inexact parameters (floats that are too small)

    Follows the approach outlined in [Puente et al. 2023](https://arxiv.org/abs/2401.00078).
"""
function isconcentrationrobust(rn::ReactionSystem; p::VarMapType = Dict())
    nps = Catalyst.get_networkproperties(rn)
    Catalyst.deficiency(rn) == 1 && !isempty(robustspecies_δ1(rn)) && return :GLOBAL_ACR

    # Check necessary condition
    possibly_ACR = false
    (length(species(rn)) == 1 || any(iszero, eachcol(complexstoichmat(rn)))) || begin
        for i in 1:length(species(rn))
            S, D = CatalystNetworkAnalysis.removespec(rn, i)
            (Catalyst.deficiency(rn) == CatalystNetworkAnalysis.deficiency(S, D)) &&
                begin
                possibly_ACR = true
                break
            end
        end
        (possibly_ACR == false) && return :NO_ACR
    end

    # Convert parameter values to rational values.
    try
        ftype, stype = eltype(p).parameters
        stype <: Rational || (p = Dict{ftype, Rational}(p))
    catch InexactError
        return :INEXACTPARAMS
    end

    # Compute steady state ideal.
    sfr_f = eval(C.SFR(rn; p = p))
    R,
        polyvars = polynomial_ring(
        QQ,
        vcat(
            map(s -> Symbolics.tosymbol(s, escape = false), species(rn)),
            map(p -> Symbolics.tosymbol(p), parameters(rn))
        )
    )
    sfr = sfr_f(polyvars...)
    specs = polyvars[1:numspecies(rn)]

    # If Grobner basis has an element of the form (x_i - α), we have ACR in i
    I = ideal(R, sfr)
    le = linearelements(I, length(specs))
    if !isempty(le)
        push!(nps.robustspecies, le...)
        return isempty(p) ? :GLOBAL_ACR : :MASS_ACTION_ACR
    end

    # Compute saturation of ideal
    J = ideal(R, prod(specs))
    satI = saturation(I, J)
    le = linearelements(satI, length(specs))
    if !isempty(le)
        push!(nps.robustspecies, le...)
        return isempty(p) ? :GLOBAL_ACR : :MASS_ACTION_ACR
    end

    # Iterate over elimination ideals I ∩ Q[x_i] and check for unique positive root.
    # Create a dummy univariate ring
    r_ξ, ξ = polynomial_ring(QQ, :ξ)

    !isempty(p) && for i in 1:numspecies(rn)
        IQ = Oscar.eliminate(I, vcat(specs[1:(i - 1)], specs[(i + 1):end]))
        for g in gens(IQ)
            iszero(g) && continue

            # Generate the univariate polynomial corresponding to the generator
            coeff_pos = [exponent_vector(g, j)[i] + 1 for j in 1:length(Oscar.terms(g))]
            coeffs = zeros(QQFieldElem, first(coeff_pos))
            coeffs[coeff_pos] .= Oscar.coefficients(g)
            poly = r_ξ(coeffs)

            iszero(poly) && continue
            if Hecke.n_positive_roots(poly) == 1
                push!(nps.robustspecies, i)
                return :MASS_ACTION_ACR
            end
        end
    end

    return :INCONCLUSIVE
end

"""
    linearelements(I, numspecies)

    Look for terms of the form x - α_i in the basis for the positive steady-state ideal. Their presence
    implies the existence of ACR in that species. 
"""
function linearelements(I::Ideal, numspecies::Int)
    G = Oscar.groebner_basis(I)
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

# TODO: Other methods.
# Compute decompositions of the positive-restriction ideal
# Numerically checking for ACR
# Compute the steady-state parameterization, check for constant components.
# Algorithm 6.1 for finding all pairs of ACR candidates

"""
    robustspecies(rn::ReactionSystem)

    For a network of deficiency one, return a vector of indices corresponding to species that are concentration robust, i.e. for every positive equilbrium, the concentration of species s will be the same. 
"""
function robustspecies_δ1(rn::ReactionSystem)
    complexes, D = reactioncomplexes(rn)
    nps = Catalyst.get_networkproperties(rn)

    if Catalyst.deficiency(rn) != 1
        error("This algorithm currently only checks for robust species in networks with deficiency one.")
    end

    # A species is concentration robust in a deficiency one network if there are two non-terminal complexes (i.e. complexes
    # belonging to a linkage class that is not terminal) that differ only in species s (i.e. their difference is some
    # multiple of s. (A + B, A) differ only in B. (A + 2B, B) differ in both A and B, since A + 2B - B = A + B).

    if !nps.checkedrobust
        tslcs = terminallinkageclasses(rn)
        Z = complexstoichmat(rn)

        # Find the complexes that do not belong to a terminal linkage class
        nonterminal_complexes = deleteat!([1:length(complexes);], vcat(tslcs...))
        robust_species = Int64[]

        for c_s in nonterminal_complexes, c_p in nonterminal_complexes

            (c_s >= c_p) && continue
            # Check the difference of all the combinations of complexes. The support is the set of indices that are non-zero
            suppcnt = 0
            supp = 0
            for i in 1:size(Z, 1)
                (Z[i, c_s] != Z[i, c_p]) && (suppcnt += 1; supp = i)
                (suppcnt > 1) && break
            end

            # If the support has length one, then they differ in only one species, and that species is concentration robust.
            (suppcnt == 1) && (supp ∉ robust_species) && push!(robust_species, supp)
        end
        nps.checkedrobust = true
        nps.robustspecies = robust_species
    end

    return nps.robustspecies
end

##### DEFICIENCY METHODS
# These are based on the computation of "robust ratios" between complexes - complexes y, y' such that x^y / x^y' is constant for all positive steady states.
#   For def 0 networks, this is true for all networks in the same linkage class
#   For def 1 networks, this is true for all
