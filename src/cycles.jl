"""
    isconsistent(rs::ReactionSystem)

    Checks if a reaction network is consistent, i.e. admits a positive equilibrium for some choice of rate constants. Equivalent to [`ispositivelydependent`](@ref).
"""
function isconsistent(rs::ReactionSystem)
    cyclemat = Catalyst.cycles(rs)
    return has_positive_solution(cyclemat)
end

"""
    isconservative(rs::ReactionSystem)

    Checks if a reaction network is conservative, i.e. admits a positive linear conserved quantity. A positive linear conserved quantity is one for which the coefficient of each species is greater than zero. 
"""
function isconservative(rs::ReactionSystem)
    conslaws = conservationlaws(rs)
    return has_positive_solution(copy(conslaws'))
end

"""
    ispositivelydependent(rs::ReactionSystem)

    See documentation for [`isconsistent`](@ref).
"""
function ispositivelydependent(rs::ReactionSystem)
    return isconsistent(rs)
end

"""
    elementary_flux_modes(rn::ReactionSystem)

    Given a reaction network, return the set of elementary flux modes of the reaction network. 
"""
function elementary_flux_modes(rn::ReactionSystem)
    S = netstoichmat(rn)
    m, n = size(S)
    hyperplanes = [Polyhedra.HyperPlane(S[i, :], 0) for i in 1:m]
    halfspaces = [Polyhedra.HalfSpace(-I(n)[i, :], 0) for i in 1:n]
    poly = Polyhedra.polyhedron(hrep(hyperplanes, halfspaces), CDDLib.Library())
    vrep(poly)

    return EFMs = Matrix{Int64}(reduce(hcat, map(x -> x.a, Polyhedra.rays(poly))))
end
