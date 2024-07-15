"""
    deficiencyonealgorithm(rn::ReactionSystem)

    Determine whether a regular deficiency one network will have the ability to admit multiple equilibria and degenerate equilibria. Returns true if so. 
"""

function deficiencyonealgorithm(rn::ReactionSystem) 
    (deficiency(rn) == 1 && isregular(rn)) || error("The deficiency one algorithm only works with regular deficiency one networks.")

    g = confluencevector(rn)
    partitions = generatepartitions(rn)
    reversible = isreversible(rn)

    # iterate over partitions
    
    for partition in partitions
        solveconstraints(g, partition) && return true
        isreversible(rn) && solveconstraints(-g, partition) && return true
    end
    
    return false
end

function isregular(rn::ReactionSystem) 
    lcs = linkageclasses(rn); tlcs = terminallinkageclasses(rn)

    # Check whether terminal linkage classes are tree-like.
    img = incidencematgraph(rn)
    tlc_graphs = [Graphs.induced_subgraph(img, tlc)[1] for tlc in tlcs] 

    isconsistent(rn) && (length(lcs) == length(tlcs)) && all(Graphs.is_tree ∘ SimpleGraph, tlc_graphs)
end

function confluencevector(rn::ReactionSystem) 
    Y = complexstoichmat(rn)
    D = incidencemat(rn)

    g = D * nullspace(Y*D)
    @assert size(g, 2) == 1

    # Check the condition for absorptive sets, round g to integer vector? 
    g = g[:, 1]
    g
end

# Given a confluence vector and a UML partition, determine whether there is a vector that is sign-compatible with the stoichiometric subspace. 
function solveconstraints(rn::ReactionSystem, confluence::Vector, partition::Array) 
    Y = complexstoichmat(rn)
    S = netstoichmat(rn)
    s, c = size(Y)
    r = size(S, 2)
    tlcs = terminallinkageclasses(rn)
    
    @variable(model, μ[1:s])
    @objective(model, Min, 0)

    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # 1) Ensure that y⋅μ = y'⋅μ for all y, y' ∈ M
    @variable(model, n)
    @constraint(model, Y[:, partition[2]]' * μ == n * ones(length(partition[2])))

    # 2) Ensure that y⋅μ > y'⋅μ whenever y is higher than y'
    @variable(model, ϵ)
    for (y_u, y_m) in IterTools.product(partition[1], partition[2])
        @constraint(model, Y[:, y_u]' * μ - Y[:, y_m]' * μ >= ϵ) 
    end

    for (y_u, y_l) in IterTools.product(partition[1], partition[3])
        @constraint(model, Y[:, y_u]' * μ - Y[:, y_l]' * μ >= ϵ) 
    end

    for (y_m, y_l) in IterTools.product(partition[2], partition[3])
        @constraint(model, Y[:, y_m]' * μ - Y[:, y_l]' * μ >= ϵ) 
    end

    # 3) Ensure cut-link condition (internal ordering within TLCs).  
    for (i, tlc) in enumerate(tlcs)
        tlc_graph, vmap = Graphs.induced_subgraph(incidencematgraph(rn), tlc)
        tlc_graph = Graphs.SimpleGraph(tlc_graph)

        for edge in edges(tlc_graph)
            sr, ds = src(edge), dst(edge)

            # For each edge, generate a 2-partition by removing the edge.
            rem_edge!(tlc_graph, sr, ds)
            ccs = Graphs.connected_components(tlc_graph)
            complexsets = [vmap[cc] for cc in ccs]

            # Identify which component the source vertex is in
            srcset = sr ∈ complexsets[1] ? 1 : 2

            # The sign of the sum of the confluence vector's components over the 
            # set of source complexes and whether y, y' lie in U or L 
            # determines whether y ⋅ μ > y' ⋅ μ or the opposite
            cutsum = sum(confluence[complexset[srcset]])
            if sr ∈ partition[1]
                @constraint(model, sign(Y[:, sr]' * μ - Y[:, ds]' * μ) == sign(cutsum)) 
            elseif src(edge) ∈ partition[3]
                @constraint(model, sign(Y[:, sr]' * μ - Y[:, ds]' * μ) == -sign(cutsum)) 
            end

            add_edge!(tlc_graph, src(edge), dst(edge))
        end
    end
    
    # 4) Check that μ is sign-compatible with the stoichiometric subspace.
    @variable(model, coeffs[1:r])
    @constraint(model, sign.(S*coeffs) == sign.(μ))

    optimize!(model)
    is_solved_and_feasible(model) 
end

function generatepartitions(rn::ReactionSystem) 
    tslcs = terminallinkageclasses(rn)
    nonterminal_complexes = deleteat!([1:numcomplexes(rn);], vcat(tslcs...))
    nt = length(tslcs)
    num_trivial_tlcs = count(tlc -> size(tlc) == 1, tlcs)

    # In this representation of the UML partition, the indexes of complexes 
    # in each component are stored. 
    partitions_complexes = Array[]

    # In this representation of the UML partition, the indexes of terminal
    # linkage classes in each component are stored. 
    partitions_tlcs = Array[]

    for i in 0:3^nt - 1
        buckets = digits(i, base=3, pad=nt)
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

        U = vcat(tslcs[partition[1]]...)
        L = vcat(tslcs[partition[2]]...)
        M = vcat(tslcs[partition[3]]..., nonterminal_complexes)
        push!(partitions_complexes, [U, L, M])
        push!(partitions_tlcs, partition)
    end
    partitions_complexes, partitions_tlcs
end
