const C = CatalystNetworkAnalysis

"""
    deficiencyonealgorithm(rn::ReactionSystem)

    Determine whether a regular deficiency one network will have the ability to admit multiple equilibria and degenerate equilibria. Returns true if so. 
"""
function deficiencyonealgorithm(rn::ReactionSystem)
    (Catalyst.deficiency(rn) == 1 && isregular(rn)) || error("The deficiency one algorithm only works for regular deficiency one networks.")

    g = confluencevector(rn)
    partitions = generatepartitions(rn)
    reversible = isreversible(rn)
    cutdict = cutlinkpartitions(rn)

    Y = complexstoichmat(rn); S = netstoichmat(rn)
    s, c = size(Y); r = size(S, 2)
    
    # Initialize sign compatibility model. 
    model = C.add_sign_constraints(S, var_name = "μ")
    @variable(model, n)

    # Iterate over partitions. For each pair of UML partition and confluence vector, 
    # check if there exists a μ that satisfies the resulting constraints. 
    for partition in partitions
        sol, feasible = solveconstraints(rn, model, g, partition, cutdict)
        feasible && return true

        if reversible
            sol, feasible = solveconstraints(rn, model, -g, partition, cutdict)
            feasible && return true
        end
    end
    
    return false
end

function isregular(rn::ReactionSystem) 
    lcs = Catalyst.linkageclasses(rn)
    tlcs = terminallinkageclasses(rn)

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

    idx = findfirst(i -> @view(g[:, i]) != zeros(nc), 1:nc) 
    g = g[:, idx]
    tlcs = terminallinkageclasses(rn)
    lcs = Catalyst.linkageclasses(rn)

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

function solveconstraints(rn::ReactionSystem, model::Model, confluence::Vector, partition::Array, cutdict::Dict; M = 1E8, ϵ = 1E-3) 
    Y = complexstoichmat(rn); S = netstoichmat(rn)
    μ = model[:μ]; n = model[:n]

    # 1) Ensure that y⋅μ = y'⋅μ for all y, y' ∈ M
    @constraint(model, M_equal, Y[:, partition[2]]' * μ == n * ones(length(partition[2])))

    # 2) Ensure that y⋅μ > y'⋅μ whenever y is higher than y'
    @constraint(model, UM[u=partition[1], m=partition[2]], Y[:, u]' * μ - Y[:, m]' * μ ≥ ϵ) 
    @constraint(model, UL[u=partition[1], l=partition[3]], Y[:, u]' * μ - Y[:, l]' * μ ≥ ϵ) 
    @constraint(model, ML[m=partition[2], l=partition[3]], Y[:, m]' * μ - Y[:, l]' * μ ≥ ϵ) 

    # 3) Ensure cut-link condition (internal ordering within TLCs).  
    for edge in keys(cutdict)
        (edge[1] ∈ partition[1] || edge[1] ∈ partition[3]) || continue
        srcset, dstset = cutdict[edge]
        sr, ds = edge
        cutsum = sum(confluence[srcset])

        if cutsum == 0
            @constraint(model, (Y[:, sr]' * μ == Y[:, ds]' * μ), base_name = "Edge $sr-$ds")
        elseif cutsum > 0
            sr ∈ partition[1] ? 
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≥ ϵ, base_name = "Edge $sr-$ds") :
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≤ -ϵ, base_name = "Edge $sr-$ds")
        elseif cutsum < 0
            sr ∈ partition[3] ? 
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≥ ϵ, base_name = "Edge $sr-$ds") : 
                @constraint(model, (Y[:, sr]' * μ - Y[:, ds]' * μ) ≤ -ϵ, base_name = "Edge $sr-$ds")
 
        end
    end

    optimize!(model) 
    feasible = is_solved_and_feasible(model)
    μ_sol = nothing 

    unregister(model, :M_equal); 
    unregister(model, :UM); unregister(model, :UL); unregister(model, :ML) 

    for con in all_constraints(model, include_variable_in_set_constraints = true)
        JuMP.name(con) == "" || begin
            unregister(model, Symbol(JuMP.name(con))) 
            delete(model, con)
        end
    end

    μ_sol, feasible
end

# Assuming that the graph is regular, generate the set of cut-link partitions
# created by removing each cut-link connecting two terminal complexes. 
function cutlinkpartitions(rn::ReactionSystem) 
    tlcs = terminallinkageclasses(rn)
    lcs = Catalyst.linkageclasses(rn)
    inclusions = [findfirst(lc -> issubset(tlc, lc), lcs) for tlc in tlcs]

    img = Graphs.SimpleGraph(incidencematgraph(rn))
    lc_graphs = [Graphs.induced_subgraph(img, lc) for lc in lcs[inclusions]]
    tlc_graphs = [Graphs.induced_subgraph(img, tlc) for tlc in tlcs]
    
    cutdict = Dict{Tuple{Int64, Int64}, Vector{Vector{Int64}}}()

    for (i, tlc) in enumerate(tlcs)
        tlc_graph, tvmap = tlc_graphs[i]
        lc_graph, lvmap = lc_graphs[i]
        tlvmap = [findfirst(==(v), lvmap) for v in tvmap] # Get the indexes of tlc in lc

        for edge in Graphs.edges(tlc_graph)
            sr, ds = Graphs.src(edge), Graphs.dst(edge)
            gsr, gds = tvmap[sr], tvmap[ds]
            lsr, lds = tlvmap[sr], tlvmap[ds]

            # For each edge, generate a 2-partition by removing the edge. TODO: probably a more efficient way to do this. 
            Graphs.rem_edge!(lc_graph, lsr, lds)
            ccs = Graphs.connected_components(lc_graph)
            complexsets = [lvmap[cc] for cc in ccs]
             
            # Identify which component the source vertex is in
            srcset, dstset = gsr ∈ complexsets[1] ? (1,2) : (2,1) 
            cutdict[(gsr, gds)] = [complexsets[srcset], complexsets[dstset]]

            Graphs.add_edge!(lc_graph, lsr, lds)
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
    numparts = 3^num_nontrivial_tlcs
    partitions_complexes = Array[]

    for i in 0:numparts - 1
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
    end
    partitions_complexes
end
