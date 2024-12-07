### Linear programming utilities for the algorithms. 

#################################
### JuMP Model Building Utils ###
#################################

"""
    add_subspace_constraints(S; model = nothing, varName = "")

    Given a matrix S, add the constraints that the solution vector μ will lie in the image space of S. If the model does not already exists, initialize one.
"""
function add_subspace_constraints(S::M; model = nothing, varName::String = "") where M <: AbstractMatrix
    (s, r) = size(S)

    # Initialize model if none provided. 
    model === nothing && begin 
        model = Model(HiGHS.Optimizer); 
        set_silent(model)
        set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)
        @objective(model, Min, 0)
    end

    coeffs = varName*"_coeffs"
    model[Symbol(coeffs)] = @variable(model, [i = 1:r], base_name = coeffs) 
    model[Symbol(varName)] = @variable(model, [i = 1:s], base_name = varName)

    @constraint(model, S*model[Symbol(coeffs)] == model[Symbol(var)])
    return model
end

"""
    add_sign_constraints(S; model = nothing, varName = "")

    Given a matrix S, add the constraints that the solution vector var will be sign-compatible with the image space of S. If the model does not already exists, initialize one.
"""
function add_sign_constraints(S::M; model = nothing, varName::String = "") where M <: AbstractMatrix
    (s, r) = size(S)

    # Initialize model if none provided. 
    model == nothing && begin 
        model = Model(HiGHS.Optimizer); 
        set_silent(model)
        set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)
        @objective(model, Min, 0)
    end
    
    ispos = var*"_ispos"; isneg = var*"_isneg"; iszer = var*"_iszero"

    model[Symbol(ispos)] = @variable(model, [i = 1:s], Bin, base_name = ispos)
    model[Symbol(isneg)] = @variable(model, [i = 1:s], Bin, base_name = isneg)
    model[Symbol(iszer)] = @variable(model, [i = 1:s], Bin, base_name = iszer) 

    @constraints(model, begin
       model[Symbol(iszer)] + model[Symbol(ispos)] + model[Symbol(isneg)] == ones(s) 
       sum(model[Symbol(iszer)]) <= s - 1 # Ensure that var is not the zero vector.

       # iszero = 1 --> var[i] == 0 <--> (S * coeffs)[i] == 1
       model[Symbol(var)] + M * (ones(s) - model[Symbol(iszer)]) ≥ zeros(s)
       model[Symbol(var)] - M * (ones(s) - model[Symbol(iszer)]) ≤ zeros(s)
       (S*model[Symbol(coeffs)]) + M * (ones(s) - model[Symbol(iszer)]) ≥ zeros(s)
       (S*model[Symbol(coeffs)]) - M * (ones(s) - model[Symbol(iszer)]) ≤ zeros(s)

       # isnegative = 1 --> var[i] < 0 <--> (S * coeffs)[i] < 0
       model[Symbol(var)] - M * (ones(s) - model[Symbol(isneg)]) ≤ -ones(s) * ϵ
       (S*model[Symbol(coeffs)]) - M * (ones(s) - model[Symbol(isneg)]) ≤ -ones(s) * ϵ

       # ispositive = 1 --> var[i] > 0 <--> (S * coeffs)[i] > 0
       model[Symbol(var)] + M * (ones(s) - model[Symbol(ispos)]) ≥ ones(s) * ϵ
       (S*model[Symbol(coeffs)]) + M * (ones(s) - model[Symbol(ispos)]) ≥ ones(s) * ϵ
    end)

    model
end

###################
### Misc. Utils ###
###################

# Given a matrix M, determine whether there is some x in the image space of M that is positive. If nonneg is true, instead looks for a non-negative solution that is not the zero vector.
function has_positive_solution(M::M{T}; nonneg = false) where {T, M <: AbstractMatrix}
    iszero(M) && return false
    all(>=(0), M) && return true
    m, n = size(M)

    model = Model(HiGHS.Optimizer)
    set_silent(model)
    @variable(model, coeffs[1:n])
    @objective(model, Min, 0)

    nonneg ? 
        @constraint(model, M * coeffs ≥ zeros(m)) : 
        @constraint(model, M * coeffs ≥ ones(m))

    optimize!(model)
    is_solved_and_feasible(model)
end

# Suppose we have a cone defined by the intersection of the nullspace of S and the positive orthant. Then is_extreme_ray tells whether a given vector in this cone is an extreme ray.
function is_extreme_ray(S::M, x::V; atol = 1e-9) where {M <: AbstractMatrix, V <: AbstractVector}
    m, n = size(S)
    length(x) != n && error("The length of x is not correct, expected $n and received $(length(x)).")
    !isapprox(S*x, zeros(m); atol) && error("The provided vector $x is not a solution to Sx = 0.")
    idxset = findall(!=(0), x)
    is_extreme_idxset(S, idxset)
end

function is_extreme_idxset(S::M, idxs::Vector{Int}) where {M <: AbstractMatrix}
    m, n = size(S)
    cone_mat = [I; S; -S] 
    cone_mat_eq  = Matrix{eltype(S)}(undef, n+2*m - length(idxs), n)

    rowidxs = deleteat!(collect(1:n+2*m), idxset)
    cone_mat_eq = @view cone_mat[rowidxs, :]
    return rank(cone_mat_eq) == n - 1
end

"""
    elementaryfluxmodes(rn::ReactionSystem)

    Given a reaction network, return the set of elementary flux modes of the reaction network. 
"""
function elementary_flux_modes(rn::ReactionSystem)
    S = netstoichmat(rn)
    m, n = size(S)
    hyperplanes = [Polyhedra.HyperPlane(S[i, :], 0) for i in 1:m]
    halfspaces = [Polyhedra.HalfSpace(-I(n)[i, :], 0) for i in 1:n]
    polycone = Polyhedra.polyhedron(hrep(hyperplanes, halfspaces))

    Polyhedra.vrep(polycone)
    Polyhedra.removevredundancy!(polycone)

    EFMs = reduce(hcat, map(x->x.a, polycone.vrep.rays.rays))
end
