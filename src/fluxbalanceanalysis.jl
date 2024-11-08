"""
    elementaryfluxmodes(rn::ReactionSystem)

    Given a reaction network, returun the set of elementary flux modes of the reaction network. 
"""
function elementaryfluxmodes(rn::ReactionSystem)
    S = netstoichmat(rn)
    m, n = size(S)
    hyperplanes = [Polyhedra.HyperPlane(S[i, :], 0) for i in 1:m]
    halfspaces = [Polyhedra.HalfSpace(-I(n)[i, :], 0) for i in 1:n]
    polycone = Polyhedra.polyhedron(hrep(hyperplanes, halfspaces))

    Polyhedra.vrep(polycone)
    Polyhedra.removevredundancy!(polycone)

    EFMs = reduce(hcat, map(x->x.a, polycone.vrep.rays.rays))
end


