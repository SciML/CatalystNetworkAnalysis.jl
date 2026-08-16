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
    l = length(Catalyst.linkageclasses(incidencematgraph(D[:, nonnull_rx])))

    return n - l - ss
end

"""
    removespec(rn, spec)

Remove one species from a reaction network's stoichiometric representation.

# Arguments

- `rn`: reaction network whose species are being reduced.
- `spec`: one-based index of the species to remove.

# Returns

A pair `(S, D)` containing the reduced net stoichiometric matrix and the
reduced incidence matrix. Complexes that become identical after removing the
species are merged before constructing `D`.
"""
function removespec(rn::ReactionSystem, spec::Int64)
    s, r = size(netstoichmat(rn))
    Y_new = @view complexstoichmat(rn)[vcat(1:(spec - 1), (spec + 1):s), :]
    complexes = transitiveclosure(eachcol(Y_new), ==)
    Y_new = @view Y_new[:, first.(complexes)]

    D_old = incidencemat(rn)
    D = SparseArrays.spzeros(Int64, length(complexes), size(D_old, 2))

    for (i, complex) in enumerate(complexes)
        row = @view D[i, :]
        row .= sum([D_old[c, :] for c in complex])
    end
    return Y_new * D, D
end

"""
    transitiveclosure(arr, relation)

Partition an iterable into connected components under a relation.

The relation is evaluated for pairs of elements and is assumed to be
symmetric and reflexive. The returned components are represented by the
one-based indices of the elements in `arr`.

# Arguments

- `arr`: iterable whose elements are to be partitioned.
- `relation`: binary predicate that identifies related elements.

# Returns

A vector of vectors of indices. Two elements occur in the same vector when
they are connected by the transitive closure of `relation`.
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
    return partition = Graphs.connected_components(intersectgraph)
end

"""
    reactiontocomplexmap(rn::ReactionSystem)

Construct a map from reaction indices to source and product complex indices.

# Arguments

- `rn`: reaction system whose incidence matrix determines the complex map.

# Returns

A `Dict{Int, Pair{Int, Int}}` mapping each reaction index to
`source_complex => product_complex`.
"""
function reactiontocomplexmap(rn::ReactionSystem)
    rxtocomplexmap = Dict{Int, Pair{Int, Int}}()
    D = incidencemat(rn)

    for i in 1:size(D, 2)
        p = findfirst(==(1), @view D[:, i])
        s = findfirst(==(-1), @view D[:, i])
        rxtocomplexmap[i] = s => p
    end
    return rxtocomplexmap
end

function matrixtree(g::SimpleDiGraph, distmx::Matrix{T}) where {T}
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

    for tree in Combinatorics.combinations(collect(edges(ug)), n - 1)
        tree = SimpleGraph(tree)
        isempty(Graphs.cycle_basis(t)) || continue

        # Add the tree product for each vertex
        for v in 1:n
            rootedTree = reverse(Graphs.bfs_tree(t, v, dir = :in))
            π[v] += treeweight(t, g, distmx)
        end
    end

    # Constructed rooted trees for every vertex, compute sum
    return π
end

function treeweight(tree::SimpleDiGraph, distmx::Matrix{T}) where {T}
    prod = 1
    for e in edges(tree)
        s = Graphs.src(e)
        p = Graphs.dst(e)
        prod *= distmx[s, p]
    end
    return prod
end

"""
    ratematrix(rs::ReactionSystem, rates = reactionrates(rs))

Construct the complex-to-complex rate matrix for a reaction system.

# Arguments

- `rs`: reaction system whose reactions define the matrix entries.
- `rates`: optional vector of reaction rate constants. When omitted, the
  symbolic reaction rates from `rs` are used.

# Returns

An `n × n` matrix, where `n` is the number of complexes and entry `(i, j)` is
the rate associated with a reaction from complex `i` to complex `j`.
"""
function ratematrix(rs::ReactionSystem, rates::Vector{T} = reactionrates(rs)) where {T}
    complexes, D = reactioncomplexes(rs)
    n = length(complexes)
    rxns = reactions(rs)
    ratematrix = (T <: SymbolicUtils.BasicSymbolic) ? zeros(Any, n, n) : zeros(T, n, n)

    for r in 1:length(rxns)
        rxn = rxns[r]
        s = findfirst(==(-1), @view D[:, r])
        p = findfirst(==(1), @view D[:, r])
        ratematrix[s, p] = rates[r]
    end
    return ratematrix
end

function matrixpower(v::Vector{T}, M::Matrix{T}) where {T <: Number}
    n, m = size(M)
    m != length(v) && error("Incorrect dimensions of matrix M.")
    out = Vector{T}(undef, n)

    for i in 1:n
        out[i] = prod([v[j]^M[i, j] for j in 1:m])
    end
    return out
end
