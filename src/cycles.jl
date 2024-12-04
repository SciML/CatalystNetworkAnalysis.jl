"""
    cycles(rs::ReactionSystem)

    Returns the matrix of cycles (or flux vectors), or reaction fluxes at steady state. These correspond to right eigenvectors of the stoichiometric matrix. Equivalent to [`fluxmodebasis`](@ref). 
"""
function cycles(rs::ReactionSystem)
    # nps = get_networkproperties(rs)
    nsm = netstoichmat(rs)
    cycles(nsm)
    # !isempty(nps.cyclemat) && return nps.cyclemat
    # nps.cyclemat = cycles(nsm; col_order = nps.col_order)
    # nps.cyclemat
end

function cycles(nsm::T; col_order = nothing) where {T<:AbstractMatrix}

    # compute the left nullspace over the integers
    N = MT.nullspace(nsm; col_order)

    # if all coefficients for a cycle are negative, make positive
    for Nrow in eachcol(N)
        all(r -> r <= 0, Nrow) && (Nrow .*= -1)
    end

    # check we haven't overflowed
    iszero(nsm * N) || error(
        "Calculation of the cycle matrix was inaccurate, " *
        "likely due to numerical overflow. Please use a larger integer " *
        "type like Int128 or BigInt for the net stoichiometry matrix.",
    )

    T(N)
end

"""
    isconsistent(rs::ReactionSystem)

    Checks if a reaction network is consistent, i.e. admits a positive equilibrium for some choice of rate constants. Equivalent to [`ispositivelydependent`](@ref).
"""
function isconsistent(rs::ReactionSystem)
    cyclemat = cycles(rs)
    has_positive_solution(cyclemat)
end

"""
    isconservative(rs::ReactionSystem)

    Checks if a reaction network is conservative, i.e. admits a positive linear conserved quantity. A positive linear conserved quantity is one for which the coefficient of each species is greater than zero. 
"""
function isconservative(rs::ReactionSystem)
    conslaws = conservationlaws(rs)
    has_positive_solution(copy(conslaws'))
end

"""
    ispositivelydependent(rs::ReactionSystem)

    See documentation for [`isconsistent`](@ref).
"""
function ispositivelydependent(rs::ReactionSystem)
    isconsistent(rs)
end

# Given a matrix M, evaluate whether there is some x in the image space of M that is positive. If nonneg is true, instead looks for a non-negative solution.  
function haspositivesolution(M::Matrix; nonneg = false) 
    isempty(M) && return false
    m, n = size(M)

    for i = 1:n
        all(>(0), @view M[:, i]) && return true
    end

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, coeffs[1:n])
    @objective(model, Min, 0)

    nonneg ? 
        @constraint(model, M * coeffs >= zeros(m)) : 
        @constraint(model, M * coeffs >= ones(m))

    optimize!(model)
    is_solved_and_feasible(model) ? true : false
end

# Suppose we have a cone defined by the intersection of the subspace Sx = 0 and the positive orthant. Then isextreme(S, x) tests whether the set of indices in x
function isextreme_poscone(S::Matrix; idxset::Vector = [], x::Vector = []) 
    m, n = size(S)
    if isempty(x) && isempty(idxset)
        error("Must provide either an index set or a solution vector.")
    elseif !isempty(x) 
        S*x == 0 && error("x is not a solution of Sx = 0.")
        idxset = findall(!=(0), x)
    end

    cone_mat = [I; S; -S] 
    cone_mat_eq = cone_mat[[deleteat!(collect(1:n), idxset)..., n+1:n+2*m...], :]
    return rank(cone_mat_eq) == n - 1
end

"""
    elementaryfluxmodes(rn::ReactionSystem)

    Given a reaction network, return the set of elementary flux modes of the reaction network. 
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
