# This file contains functionality for translated reaction networks, which extend the results
# of deficiency theory to a much wider range of networks. 

# (1) Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.
 
# Functionality for translating chemical reaction networks. 

"""
    translate(rn::ReactionNetwork)

    Given a reaction network, attempt to find a strongly-resolvable translation that is weakly-reversible and deficiency zero. Such translations are useful because they can be easily parameterized, and their multistability characteristics can be easily determined. Follows the algorithm implemented in (Johnston, 2018). 
"""
# Eventually we want this to just construct the lowest-deficiency network, see Johnston/Tonello 2017
function WRDZ_translation(rn::ReactionSystem) 
    Y_K = complexstoichmat(rn)
    Y_T = zeros(size(kineticcomplexes))
    e = collect(Graphs.edges(incidencematgraph(rn)))

    rr_adj = rrgraph(rn)
    r = length(reactions(rn))
    
    P1 = collect(1:r); P2 = Vector{Int64}(); P3 = Vector{Int64}()
    while true
        @label step_2
        i = first(P1)
        P2 = [r]; popfirst!(P1)

        # Update reactant complexes
        @label step_3
        rxs = findall(==(1), rr_adj[i, :])
        for j in rxs
            p, s = (dst(e[i]), src(e[j]))
            (@view Y_T[:, j]) .= Y_K[:, p] - Y_K[:, s] + Y_T[:, i]
            P3 = setdiff(push!(P3, j), P2)
        end

        !isempty(P3) && begin
            P2 = copy(P3); P1 = setdiff(P1, P2); empty!(P3)
            @goto step_3
        end

        !isempty(P1) && @goto step_2
        isempty(P3) && isempty(P1) && break
    end

    Y_T
end

function rrgraph(rn::ReactionSystem) 
    satisfiesdeficiencyzero(rn) && return nothing

    # Get the complexes and define generalized complexes
    kineticcomplexes = complexstoichmat(rn)
    translatedcomplexes = copy(complexstoichmat(rn))

    # Find elementary flux modes and reactions with shared source complexes
    EFMs = elementaryfluxmodes(rn)
    # Check that EFMs are unitary and cover R

    # Solve BLP
    model = Model(HiGHS.Optimizer);
    set_silent(model)
    set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)

    # Essentially an adjacency matrix for the reaciton-reaction graph 
    @variable(model, edge[1:r, 1:r], Bin)
    @objective(model, Min, sum(edge))

    # Add CS, EM, PS constraints
    CS, EM = rrgraphpartitions(rn, EFMs)
    # Partition constraint
    for i in 1:length(EM), j in 1:length(EM)
        (i == j || i > j) && continue
        @constraint(model, edge[EM[i], EM[j]] .== 0)
        @constraint(model, edge[EM[j], EM[i]] .== 0)
    end

    # CS constraints: any two reactions with the same source complex have the same incoming edges 
    for part in CS 
        for i in 1:length(part), j in 1:length(part)
            (i == j || i > j) && continue
            @constraint(model, edge[:, part[i]] - edge[:, part[j]] .== 0)
        end
    end

    # EM constraints
    supp = Vector{Int64}()
    for col in eachcol(EFMs)
        supp = findall(>(0), col)
        l = length(supp)
        @constraint(model, sum(edge[supp, supp]) == l)
        for i in 1:length(supp), j in 1:length(supp)
            (i == j || i > j) && continue
            @constraint(model, sum(edge[:, j]) == 1)
            @constraint(model, sum(edge[i, :]) == 1)
        end

        # Prevent subcycles. 
        for i in 2:floor(Int, l / 2)
            for idxs in combinations(supp, i)
                @constraint(model, sum(edge[idxs, idxs]) <= i - 1)
            end
        end
    end

    optimize!(model)
    is_solved_and_feasible(model) && return edge
end

# Generate partitions of the reaction-to-reaction graph
function rrgraphpartitions(rn::ReactionSystem, efms) 
    D = incidencemat(rn)
    CSrxs = Vector{Vector{Int64}}()
    outrxs = Vector{Int64}()
    for i in 1:size(D, 1)
        outrxs = findall(==(-1), @view D[i, :])
        length(outrxs) > 1 ? push!(CSrxs, outrxs) : continue
    end

    EFMsupports = deepcopy(CSrxs)
    EFMpartitions = [Vector{Int64}() for _ in 1:length(EFMsupports)]
    for j in 1:size(efms, 2)
        supprxs = findall(>(0), @view efms[:, j])
        for (i, supp) in enumerate(EFMsupports)
            any(∈(supp), supprxs) && begin
                push!(EFMpartitions[i], j)
                for v in supprxs
                    (v ∉ supp) && push!(supp, v)
                end
                continue
            end
        end
    end
    # Collapse common elements of the partition
    ccs = transitiveclosure_intersect(EFMsupports)
    EFMsupports = [union(EFMsupports[cc]...) for cc in ccs]
    EFMpartitions = [union(EFMpartitions[cc]...) for cc in ccs]
    CSrxs, EFMpartitions, EFMsupports
end

# Compute the vectors that form a component under transitive closure under intersection
function transitiveclosure_intersect(vectors::Vector{T}) where {T <: Union{Vector, Set}}
    adjmat = spzeros(Bool, length(vectors), length(vectors))
    for idx in CartesianIndices(adjmat)
        (i, j) = Tuple(idx)
        if i == j
            continue
        elseif i > j
            adjmat[idx] = !isdisjoint(vectors[i], vectors[j])
        else
            adjmat[idx] = adjmat[j, i]
        end
    end

    intersectgraph = SimpleGraph(adjmat)
    cc = Graphs.connected_components(intersectgraph)
end


# (1) Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).

# Generating the parameterization from the above paper. 

# Kinetic deficiency zero reaction networks have positive parameterizations. 
function kineticdeficiency(rn::ReactionSystem) 
end

"""
    Given a reaction system, compute the positive parameterization of the system. 
    1) For generalized mass action systems with kinetic deficiency zero. 
    2) Algorithm for positive kinetic deficiency systems. 
    3) Monomial parameterization for systems with toric steady states (Millan et al. 2012)
    4) Using Matroids (fall-back), checking subsets of d variables? 

    The output of this function is a symbolic function representing the parameterization. It takes some subset of the species and maps them into a steady state.  
"""
function positiveparameterization(rn::ReactionSystem; variable_search = false) 
    if kineticdeficiency(rn) == 0
    end

    # κ ∘ τ
    # Find a spanning forest
    # Find B such that im B = ker M^T
    # H is a generalized inverse of M
end

