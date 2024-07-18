"""
    deficiencyonealgorithm(rn::ReactionSystem)

    Determine whether a regular deficiency one network will have the ability to admit multiple equilibria and degenerate equilibria. Returns true if so. 
"""

function deficiencyonealgorithm(rn::ReactionSystem) 
    (deficiency(rn) == 1 && isregular(rn)) || error("The deficiency one algorithm only works with regular deficiency one networks.")

    g = confluencevector(rn)
    partitions = generatepartitions(rn)
    reversible = isreversible(rn)

    # Iterate over partitions. For each pair of UML partition and confluence vector, 
    # check if there exists a μ that satisfies the resulting constraints. 
    
    for partition in partitions
        solveconstraints(rn, g, partition) && return true
        isreversible(rn) && solveconstraints(rn, -g, partition) && return true
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
    D = incidencemat(rn)

    g = D * nullspace(Y*D)

    # Check the condition for absorptive sets, round g to integer vector? 
    g = g[:, 1]
    g = [isapprox(0, g_i, atol=1e-12) ? 0 : g_i for g_i in g]
    # min = minimum(filter(!=(0), g))
    # g ./= min
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
function solveconstraints(rn::ReactionSystem, confluence::Vector, partition::Array) 
    s, c = size(complexstoichmat(rn)); r = size(netstoichmat(rn), 2)
    
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    @variable(model, μ[1:s])
    @objective(model, Min, 0)

    # 1) Ensure that y⋅μ = y'⋅μ for all y, y' ∈ M
    @variable(model, n)
    @constraint(model, Y[:, partition[2]]' * μ == n * ones(length(partition[2])))

    # 2) Ensure that y⋅μ > y'⋅μ whenever y is higher than y'
    for (y_u, y_m) in IterTools.product(partition[1], partition[2])
        @constraint(model, Y[:, y_u]' * μ - Y[:, y_m]' * μ >= eps()) 
    end

    for (y_u, y_l) in IterTools.product(partition[1], partition[3])
        @constraint(model, Y[:, y_u]' * μ - Y[:, y_l]' * μ >= eps()) 
    end

    for (y_m, y_l) in IterTools.product(partition[2], partition[3])
        @constraint(model, Y[:, y_m]' * μ - Y[:, y_l]' * μ >= eps()) 
    end

    # 3) Ensure cut-link condition (internal ordering within TLCs).  
    tlcs = terminallinkageclasses(rn)
    for (i, tlc) in enumerate(tlcs)
        tlc_graph, vmap = Graphs.induced_subgraph(incidencematgraph(rn), tlc)
        tlc_graph = Graphs.SimpleGraph(tlc_graph)

        for edge in Graphs.edges(tlc_graph)
            sr, ds = Graphs.src(edge), Graphs.dst(edge)

            # For each edge, generate a 2-partition by removing the edge.
            Graphs.rem_edge!(tlc_graph, sr, ds)
            ccs = Graphs.connected_components(tlc_graph)
            complexsets = [vmap[cc] for cc in ccs]

            # Identify which component the source vertex is in
            srcset = sr ∈ complexsets[1] ? 1 : 2

            # The sign of the sum of the confluence vector's components over the 
            # set of source complexes and whether y, y' lie in U or L 
            # determines whether y⋅μ > y'⋅μ or the opposite
            cutsum = sum(confluence[complexsets[srcset]])
            if sr ∈ partition[1]
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ)*(cutsum) >= ϵ) 
            elseif src(edge) ∈ partition[3]
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ)*(cutsum) <= ϵ) 
            end

            Graphs.add_edge!(tlc_graph, Graphs.src(edge), Graphs.dst(edge))
        end
    end
    
    # 4) Check that μ is sign-compatible with the stoichiometric subspace.
    @variable(model, coeffs[1:r])
    @variable(model, ispositive[1:s], Bin)
    @variable(model, isnegative[1:s], Bin)
    @variable(model, iszero[1:s], Bin)
    @variable(model, isposprod[1:s], Bin)
    const M = 2^16 

    ### The logical constraints are: 
    #   iszero <--> ¬isposprod
    #   ispositive --> isposprod
    #   isnegative --> isposprod 
    #   Exactly one of ispositive, iszero, isnegative
    #
    #   if μ[i] is not zero, then it will be positive or negative, meaning that the
    #   corresponding vector in the stoichiometric subspace will be positive or
    #   negative, meaning that their product at that index will be positive. 
    @constraint(model, iszero + isposprod == ones(s))
    @constraint(model, iszero + ispositive + isnegative == ones(s))
    @constraint(model, isposprod - isnegative >= zeros(s))
    @constraint(model, isposprod - ispositive >= zeros(s))

    #   iszero == 1 --> μ[i] = 0
    @constraint(model, μ + M * (ones(s) - iszero) ≥ zeros(s))
    @constraint(model, μ - M * (ones(s) - iszero) ≤ zeros(s))
    #   isnegative == 1 --> μ[i] < 0
    @constraint(model, μ - M * (ones(s) - isnegative) ≤ ones(s) * eps())
    #   ispositive == 1 --> μ[i] > 0
    @constraint(model, μ + M * (ones(s) - ispositive) ≥ ones(s) * eps())

    #   isposprod == 1 --> (S * coeffs)[i] * μ[i] > 0
    @constraint(model, (S*coeffs) .* μ ≥ isposprod * eps())
    #   isposprod == 0 --> (S * coeffs)[i] * μ[i] <= 0
    @constraint(model, (S*coeffs) .* μ - M * isposprod ≤ 0)

    optimize!(model)
    is_solved_and_feasible(model) 
end

function generatepartitions(rn::ReactionSystem) 
    tslcs = terminallinkageclasses(rn)
    numcomplexes = size(incidencemat(rn))[1]
    nonterminal_complexes = deleteat!([1:numcomplexes;], vcat(tslcs...))
    nt = length(tslcs)
    num_trivial_tlcs = count(tlc -> length(tlc) == 1, tslcs)

    # In this representation of the UML partition, the indexes of complexes 
    # in each component are stored. 
    partitions_complexes = Array[]

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
        # push!(partitions_tlcs, partition)
    end
    partitions_complexes
end
