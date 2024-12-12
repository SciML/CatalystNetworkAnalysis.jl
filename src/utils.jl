function deficiency(S::Matrix{Int64}, D::Matrix{Int64})
    n = size(S, 1)
    ss = rank(S)
    l = length(linkageclasses(incidencematgraph(D)))

    return n - l - ss 
end

function deficiency(S::SparseMatrixCSC{Int, Int}, D::SparseMatrixCSC{Int, Int})
    n = size(D, 1)
    ss = rank(Matrix(S))
    nonnull_rx = findall(!iszero, eachcol(D))
    l = length(linkageclasses(incidencematgraph(D[:, nonnull_rx])))

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

"""
    function reactiontocomplexmap(rn::ReactionSystem)

Construct a map from the reactions of the system to the indices of the complexes they transform between.
"""
function reactiontocomplexmap(rn::ReactionSystem) 
    rxtocomplexmap = Dict{Int, Pair{Int, Int}}()
    D = incidencemat(rn)

    for i in 1:size(D, 2)
        p = findfirst(==(1), @view D[:, i])
        s = findfirst(==(-1), @view D[:, i])
        rxtocomplexmap[i] = s => p
    end
    rxtocomplexmap
end

function matrixtree(g::SimpleDiGraph, distmx::Matrix)
    n = nv(g)
    if size(distmx) != (n, n)
        error("Size of distance matrix is incorrect.")
    end

    π = zeros(n)

    if !Graphs.is_connected(g)
        ccs = Graphs.connected_components(g)
        for cc in ccs
            sg, vmap = Graphs.induced_subgraph(g, cc)
            distmx_s = distmx[cc, cc]
            π_j = matrixtree(sg, distmx_s)
            π[cc] = π_j
        end
        return π
    end

    # generate all spanning trees
    ug = SimpleGraph(SimpleDiGraph(g))
    trees = collect(Combinatorics.combinations(collect(edges(ug)), n - 1))
    trees = SimpleGraph.(trees)
    trees = filter!(t -> isempty(Graphs.cycle_basis(t)), trees)

    # constructed rooted trees for every vertex, compute sum
    for v in 1:n
        rootedTrees = [reverse(Graphs.bfs_tree(t, v, dir = :in)) for t in trees]
        π[v] = sum([treeweight(t, g, distmx) for t in rootedTrees])
    end

    # sum the contributions
    return π
end

function treeweight(t::SimpleDiGraph, g::SimpleDiGraph, distmx::Matrix)
    prod = 1
    for e in edges(t)
        s = Graphs.src(e)
        t = Graphs.dst(e)
        prod *= distmx[s, t]
    end
    prod
end
