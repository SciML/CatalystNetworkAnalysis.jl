### Linear programming utilities. 

# Given a matrix M, determine whether there is some x in the image space of M that is positive. If nonneg is true, instead looks for a non-negative solution.  
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
