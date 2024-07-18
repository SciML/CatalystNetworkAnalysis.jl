"""
    deficiencyonealgorithm(rn::ReactionSystem)

    Determine whether a regular deficiency one network will have the ability to admit multiple equilibria and degenerate equilibria. Returns true if so. 
"""

function deficiencyonealgorithm(rn::ReactionSystem) 
    (deficiency(rn) == 1 && isregular(rn)) || error("The deficiency one algorithm only works with regular deficiency one networks.")

    g = confluencevector(rn)
    partitions = generatepartitions(rn)
    reversible = isreversible(rn)
    cutidct = cutlinkpartitions(rn)

    # Iterate over partitions. For each pair of UML partition and confluence vector, 
    # check if there exists a μ that satisfies the resulting constraints. 
    
    Y = complexstoichmat(rn); S = netstoichmat(rn)
    s, c = size(Y); r = size(S, 2)
    
    # Initialization
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    M = 1E8; ϵ = 1E-8
    @variable(model, μ[1:s])
    @objective(model, Min, 0)

    for partition in partitions
        solveconstraints(rn, g, partition, cutdict)[2] && return true
        isreversible(rn) && solveconstraints(rn, -g, partition, cutdict)[2] && return true
    end
    
    return false
end

function isregular(rn::ReactionSystem) 
    lcs = linkageclasses(rn); tlcs = terminallinkageclasses(rn)

    img = incidencematgraph(rn)
    tlc_graphs = [Graphs.induced_subgraph(img, tlc)[1] for tlc in tlcs] 

    # The three conditions for regularity:
    #   1) The reaction network is consistente
    #   2) The reaction network has one terminal LC for each LC
    #   3) Each TLC is tree-like, i.e. every reaction link y <--> y' is a cut-link
    isconsistent(rn) && (length(lcs) == length(tlcs)) && all(Graphs.is_tree ∘ SimpleGraph, tlc_graphs)
end

function confluencevector(rn::ReactionSystem) 
    Y = complexstoichmat(rn)
    D = incidencemat(rn); nc = size(Y)[2]

    numcols, kerYD = nullspace_right_rational(ZZMatrix(Y*D))
    g = D*Matrix{Int}(kerYD)
    @assert rank(g) == 1

    # Check the condition for absorptive sets, round g to integer vector? 
    idx = findfirst(i -> @view(g[:, i]) != zeros(nc), 1:nc) 
    g = g[:, idx]
    tlcs = terminallinkageclasses(rn); lcs = linkageclasses(rn)

    absorptive = true
    for tlc in tlcs
        is_union = any(==(tlc), lcs) 
        if !is_union && sum(g[tlc]) < 0
            absorptive = false
            break
        end
    end

    absorptive ? g : -g
end

# Given a confluence vector and a UML partition, determine whether there is a vector that is sign-compatible with the stoichiometric subspace. 
function solveconstraints(rn::ReactionSystem, confluence::Vector, partition::Array, cutdict::Dict) 
    Y = complexstoichmat(rn); S = netstoichmat(rn)
    s, c = size(Y); r = size(S, 2)
    
    # Initialization
    model = Model(HiGHS.Optimizer)
    set_silent(model)
    M = 1E8; ϵ = 1E-8
    @variable(model, μ[1:s])
    @objective(model, Min, 0)

    # 1) Ensure that y⋅μ = y'⋅μ for all y, y' ∈ M
    @variable(model, n)
    @constraint(model, M_equal, Y[:, partition[2]]' * μ == n * ones(length(partition[2])))

    # 2) Ensure that y⋅μ > y'⋅μ whenever y is higher than y'
    @constraint(model, UM[u=partition[1], m=partition[2]], Y[:, u]' * μ - Y[:, m]' * μ ≥ ϵ) 
    @constraint(model, ML[m=partition[2], l=partition[3]], Y[:, m]' * μ - Y[:, l]' * μ ≥ ϵ) 

    # 3) Ensure cut-link condition (internal ordering within TLCs).  
    for edge in keys(cutdict)
        (edge[1] ∈ partition[1] || edge[1] ∈ partition[3]) || continue
        srcset, dstset = cutdict[edge]
        sr, ds = edge
        cutsum = sum(confluence[srcset])

        if cutsum == 0
            @constraint(model, (Y[:, sr]' * μ == Y[:, ds]' * μ))
        elseif cutsum > 0
            sr ∈ partition[1] ? 
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≥ ϵ) :
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≤ -ϵ) 
        elseif cutsum < 0
            sr ∈ partition[3] ? 
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≥ ϵ) :
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≤ -ϵ) 
        end
    end

    # 4) Check that μ is sign-compatible with the stoichiometric subspace.
    @variable(model, coeffs[1:r])
    @variable(model, ispositive[1:s], Bin)
    @variable(model, isnegative[1:s], Bin)
    @variable(model, iszero[1:s], Bin)
    M = 2^16 

    ### The logical constraints are: 
    #   Exactly one of ispositive, iszero, isnegative
    @constraint(model, iszero + ispositive + isnegative == ones(s))

    #   iszero == 1 --> μ[i] = 0 <--> S*coeffs[i] = 0
    @constraint(model, μ + M * (ones(s) - iszero) ≥ zeros(s))
    @constraint(model, μ - M * (ones(s) - iszero) ≤ zeros(s))
    @constraint(model, (S*coeffs) + M * (ones(s) - iszero) ≥ zeros(s))
    @constraint(model, (S*coeffs) - M * (ones(s) - iszero) ≤ zeros(s))
    #   isnegative == 1 --> μ[i] < 0 <--> S*coeffs[i] < 0
    @constraint(model, μ - M * (ones(s) - isnegative) ≤ ones(s) * ϵ)
    @constraint(model, (S*coeffs) - M * (ones(s) - isnegative) ≤ ones(s) * ϵ)
    #   ispositive == 1 --> μ[i] > 0 <--> S*coeffs[i] > 0
    @constraint(model, μ + M * (ones(s) - ispositive) ≥ ones(s) * ϵ)
    @constraint(model, (S*coeffs) + M * (ones(s) - ispositive) ≥ ones(s) * ϵ)

    optimize!(model)
    JuMP.value.(μ), is_solved_and_feasible(model) 
end

function cutlinkpartitions(rn::ReactionSystem) 
    tlcs = terminallinkageclasses(rn)
    cutdict = Dict{Tuple{Int64, Int64}, Vector{Vector{Int64}}}()
    
    for (i, tlc) in enumerate(tlcs)
        tlc_graph, vmap = Graphs.induced_subgraph(incidencematgraph(rn), tlc)
        tlc_graph = Graphs.SimpleGraph(tlc_graph)

        for edge in Graphs.edges(tlc_graph)
            # For each edge, generate a 2-partition by removing the edge. TODO: probably a more efficient way to do this. 
            Graphs.rem_edge!(tlc_graph, edge)
            ccs = Graphs.connected_components(tlc_graph)
            complexsets = [vmap[cc] for cc in ccs]
             
            # Identify which component the source vertex is in
            sr, ds = vmap[Graphs.src(edge)], vmap[Graphs.dst(edge)]
            srcset, dstset = sr ∈ complexsets[1] ? (1,2) : (2,1) 
            cutdict[(sr, ds)] = [complexsets[srcset], complexsets[dstset]]

            Graphs.add_edge!(tlc_graph, Graphs.src(edge), Graphs.dst(edge))
        end
    end
    cutdict
end

# Generate partition of the reactive complexes. 
function generatepartitions(rn::ReactionSystem) 
    tlcs = terminallinkageclasses(rn)
    nontrivial_tlcs = filter(tlc -> length(tlc) != 1, tlcs)
    num_nontrivial_tlcs = length(nontrivial_tlcs)
    num_trivial_tlcs = length(tlcs) - num_nontrivial_tlcs 

    num_complexes = size(incidencemat(rn))[1]
    nonterminal_complexes = deleteat!([1:num_complexes;], vcat(tlcs...))

    # In this representation of the UML partition, the indexes of complexes 
    # in each component are stored. 
    partitions_complexes = Array[]

    for i in 0:3^num_nontrivial_tlcs - 1
        buckets = digits(i, base=3, pad=num_nontrivial_tlcs)
        partition = [findall(==(l), buckets) for l in 0:2] 

        # We only keep partitions where U is lexicographically less than L,
        # since inversions where U = L', L = U' don't need to be checked twice. 
        (partition[1] <= partition[3]) || continue

        # If the number of trivial terminal linkage classes is zero, 
        # we can ignore partitions that have empty U or L
        num_trivial_tlcs == 0 && (isempty(partition[1]) || isempty(partition[3])) && continue

        # If the number of trivial terminal linkage classes is one, we can ignore
        # the partition for which both U and L are empty
        num_trivial_tlcs == 1 && (isempty(partition[1]) && isempty(partition[3])) && continue

        U = reduce(vcat, nontrivial_tlcs[partition[1]], init = Int64[])
        M = reduce(vcat, [nontrivial_tlcs[partition[2]]..., nonterminal_complexes], init = Int64[])
        L = reduce(vcat, nontrivial_tlcs[partition[3]], init = Int64[])
        push!(partitions_complexes, [U, M, L])
        # push!(partitions_tlcs, partition)
    end
    partitions_complexes
end
