function deficiency(S::Matrix{Int64}, D::Matrix{Int64})
    n = size(S, 1)
    ss = rank(S)
    l = length(linkageclasses(incidencematgraph(D)))

    return n - l - ss 
end

function deficiency(S::SparseMatrixCSC{Int, Int}, D::SparseMatrixCSC{Int, Int})
    n = size(D, 1)
    ss = rank(Matrix(S))
    l = length(linkageclasses(incidencematgraph(D)))

    return n - l - ss 
end

"""
    removespec(rn, spec)

    Return the net stoichiometric matrix and incidence matrix obtained when removing a species at index spec from the reaction network. 
"""
function removespec(rn::ReactionSystem, spec::Int64) 
    s, r = size(netstoichmat(rn))
    Y_new = @view complexstoichmat(rn)[vcat(1:spec-1, spec+1:s), :]
    complexes = transitiveclosure(eachcol(Y_new), ==)
    Y_new = @view Y_new[:, first.(complexes)]

    D_old = incidencemat(rn)
    D = SparseArrays.spzeros(Int64, length(complexes), size(D_old, 2))

    for (i, complex) in enumerate(complexes)
        row = @view D[i, :]
        row .= sum([D_old[c, :] for c in complex])
    end
    Y_new * D, D
end

"""
    transitiveclosure(arr, relation)

    Given an iterable array and a (symmetric, reflexive) relation (a function that takes two objects, and returns true if they belong to the relation), return the partitions of the set under the transitive closure of the given relation. 
"""

function transitiveclosure(arr, relation) 
    adjmat = SparseArrays.spzeros(Bool, length(arr), length(arr))
    for idx in CartesianIndices(adjmat)
        (i, j) = Tuple(idx)
        if i == j
            continue
        elseif i > j
            adjmat[idx] = relation(arr[i], arr[j])
        else
            adjmat[idx] = adjmat[j, i]
        end
    end

    intersectgraph = SimpleGraph(adjmat)
    partition = Graphs.connected_components(intersectgraph)
end
