"""
    isconsistent(rs::ReactionSystem)

Check whether a reaction system is consistent.

A consistent reaction system admits a positive equilibrium for some choice of
rate constants.

# Arguments

- `rs`: reaction system to analyze.

# Returns

`true` when the stoichiometric cycle matrix admits a positive solution and
`false` otherwise.
"""
function isconsistent(rs::ReactionSystem)
    cyclemat = Catalyst.cycles(rs)
    return has_positive_solution(cyclemat)
end

"""
    isconservative(rs::ReactionSystem)

Check whether a reaction system is conservative.

A conservative reaction system admits a positive linear conserved quantity,
meaning that every species has a strictly positive coefficient.

# Arguments

- `rs`: reaction system to analyze.

# Returns

`true` when a positive conservation law exists and `false` otherwise.
"""
function isconservative(rs::ReactionSystem)
    conslaws = conservationlaws(rs)
    return has_positive_solution(copy(conslaws'))
end

"""
    ispositivelydependent(rs::ReactionSystem)

Alias for [`isconsistent`](@ref). A reaction system is positively dependent
when its stoichiometric cycle matrix admits a positive solution.

# Arguments

- `rs`: reaction system to analyze.

# Returns

A `Bool` indicating whether the system is positively dependent.
"""
function ispositivelydependent(rs::ReactionSystem)
    return isconsistent(rs)
end

"""
    elementary_flux_modes(rn::ReactionSystem)

Compute the elementary flux modes of a reaction system.

# Arguments

- `rn`: reaction system whose stoichiometric cone is analyzed.

# Returns

An integer matrix whose columns are the elementary flux modes.
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
