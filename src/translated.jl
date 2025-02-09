# This file contains functionality for translated reaction networks, which extend the results
# of deficiency theory to a much wider range of networks. 

# (1) Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.
 
# Struct representing a network translation.
mutable struct Translation{T <: Integer}
    """The original reaction network."""
    rn::ReactionSystem
    """The complex-composition matrix of the translated reaction network."""
    Y_T::Matrix{T} 
    """The incidence matrix of the translated reaction network."""
    D_T::Matrix{T}
    """Vector X where X[i] gives the index of the complex that i was translated into."""
    translatedcomplexmap::Vector{T}
    checkedrev::Bool = false
    isweaklyrev::Bool = true
    linkageclasses::Vector{Vector{T}} = Vector{Vector{T}}(undef, 0)
    effectivedeficiency::T = -1
    kineticdeficiency::T = -1
end

"""
    WRDZ_translation(rn::ReactionNetwork)

    Given a reaction network, attempt to find a strongly-resolvable translation that is weakly-reversible and deficiency zero. Such translations are useful because they can be easily parameterized, and their multistability characteristics can be easily determined. Follows the algorithm implemented in (Johnston, 2018). 

    Returns `nothing` if no translation is possible.
"""
function WRDZ_translation(rn::ReactionSystem) 
    rr_adj = construct_rr_graph(rn)
    isnothing(rr_adj) && return nothing

    Y_K = complexstoichmat(rn)
    Y_T = copy(Y_K)
    D = incidencemat(rn, sparse = true)

    nums = size(Y_K, 1)
    (numc, numr) = size(D)
    Λ = zeros(Int, nums, numr) # the set of translation complexes
    
    P1 = Set(collect(1:numr))
    P2 = Set{Int64}()
    P3 = Set{Int64}()

    while true
        @label step_2
        i = pop!(P1)
        empty!(P2)
        push!(P2, i)

        # Update reactant complexes
        @label step_3
        rxs = findall(==(1), @view rr_adj[i, :])
        for j in rxs
            p = findfirst(==(1), @view D[:, i])
            s = findfirst(==(-1), @view D[:, j])

            for k in 1:size(Y_T, 1)
                Λ[k, j] = Y_K[k, p] - Y_K[k, s] + Λ[k, i]
            end
            setdiff!(union!(P3, j), P2)
        end

        !isempty(P3) && begin
            P2 = copy(P3)
            setdiff!(P1, P2)
            empty!(P3)
            @goto step_3
        end

        !isempty(P1) && @goto step_2
        isempty(P3) && isempty(P1) && break
    end

    rcmap = CatalystNetworkAnalysis.reactiontocomplexmap(rn)
    for i in 1:numr
        s = rcmap[i].first
        p = rcmap[i].second

        @. Y_T[:,s] += Λ[:,i]
        @. Y_T[:,p] += Λ[:,i]
    end

    # If any translated complex has negative amounts of some species, add a complex to each member of its linkage class 
    any(<(0), Y_T) && begin
        lcs = Catalyst.linkageclasses(rn)

        for lc in lcs
            complexes = @view Y_T[:, lc]
            negidxs = findall(<(0), complexes)
            negvals = complexes[negidxs]
            specidxs = map(x->x.I[1], negidxs)

            addcomplex = zeros(Int, nums)
            addcomplex[specidxs] .= -negvals

            for c in lc
                @. Y_T[:,c] += addcomplex
            end
        end
    end

    _Y_T = reduce(hcat, unique(eachcol(Y_T)))
    translatedcmap = zeros(Int, numc)
    for i in 1:numc
        j = findfirst(==(Y_T[:, i]), eachcol(_Y_T))
        translatedcmap[i] = j 
    end

    D_T = zeros(Int, size(_Y_T, 2), numr)
    for i in 1:numr
        s = rcmap[i].first
        p = rcmap[i].second
        D_T[translatedcmap[s], i] = -1
        D_T[translatedcmap[p], i] = 1
    end

    return Translation(rn = rn, Y_T = _Y_T, D_T = D_T, translatedcomplexmap = translatedcmap)
end

# Construct a reaction-reaction graph that is common-source compatible and elementary mode compatible. Return an adjacency matrix.
function construct_rr_graph(rn::ReactionSystem) 
    # Do nothing for already WRDZ networks
    Catalyst.satisfiesdeficiencyzero(rn) && return nothing

    # Find elementary flux modes and reactions with shared source complexes
    EFMs = elementary_flux_modes(rn)

    # Check that EFMs are unitary and cover R
    all(i -> (i==1 || i==0), EFMs) || return nothing
    all(>(0), sum(EFMs, dims=2)) || return nothing

    EFM_supports = [findall(>(0), efm) for efm in eachcol(EFMs)]

    model = Model(HiGHS.Optimizer);
    set_silent(model)
    set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)

    r = length(reactions(rn))
    @variable(model, edge[1:r, 1:r], Bin)
    @objective(model, Min, sum(edge))

    for i in 1:r
        @constraint(model, edge[i, i] == 0)
    end

    # Compute common sources and partitions
    csrs = common_source_reactions(rn)
    efm_parts = efm_partitions(rn, EFM_supports, csrs)

    # Partition constraint: no edges between different partitions.
    for i in 1:length(efm_parts), j in 1:length(efm_parts)
        (i == j || i > j) && continue
        @constraint(model, edge[efm_parts[i], efm_parts[j]] .== 0)
        @constraint(model, edge[efm_parts[j], efm_parts[i]] .== 0)
    end

    # CS constraints: any two reactions with the same source complex have the 
    # same incoming edges in the RR graph.
    for rxs in csrs
        for i in 1:length(rxs), j in 1:length(rxs)
            (i == j || i > j) && continue
            @constraint(model, edge[:, rxs[i]] - edge[:, rxs[j]] .== 0)
        end
    end

    # EM constraints: every EFM corresponds to a directed cycle in the RR graph.
    for (i, supp) in enumerate(EFM_supports)
        l = length(supp)

        @constraint(model, sum(edge[supp, supp]) == l)
        for r in supp
            @constraint(model, sum(edge[supp, r]) == 1)
            @constraint(model, sum(edge[r, supp]) == 1)
        end

        # Prevent subcycles. 
        for n in 2:floor(Int, l / 2)
            for idxs in combinations(supp, n)
                @constraint(model, sum(edge[idxs, idxs]) <= n - 1)
            end
        end
    end

    optimize!(model)
    is_solved_and_feasible(model) ? (return Matrix{Int64}(JuMP.value.(model[:edge]))) : (return nothing)
end

function common_source_reactions(rn::ReactionSystem) 
    D = incidencemat(rn)
    common_source_rxs = Vector{Vector{Int64}}()
    outgoing_rxs = Vector{Int64}()

    for i in 1:size(D, 1)
        outgoing_rxs = findall(==(-1), @view D[i, :])
        length(outgoing_rxs) > 1 ? push!(common_source_rxs, outgoing_rxs) : continue
    end
    common_source_rxs
end

# 1) any two EFMs with a reaction in common are in the same partition, and 2) any two EFMs with reactions with the same source complex are in the same partition.
function is_equiv_efm(efm1, efm2, csrs)
    !isdisjoint(efm1, efm2) && return true

    for rs in csrs
        (isdisjoint(efm1, rs) || isdisjoint(efm2, rs)) && continue
        !isdisjoint(efm1 ∪ rs, efm2) && return true
    end

    return false
end

# Generate partitions of the reactions such that any two reactions in the same class are in the same EFM partition
function efm_partitions(rn::ReactionSystem, efm_supports, csrs) 
    # Collapse common elements of the partition

    ccs = CatalystNetworkAnalysis.transitiveclosure(efm_supports, (efm1, efm2) -> is_equiv_efm(efm1, efm2, csrs))
    efm_parts = [union(efm_supports[cc]...) for cc in ccs]
end

# (1) Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).

# Network analysis for translations

function isweaklyreversible(trn::Translation) 
    trn.checkedrev && return trn.isweaklyrev

    g = Catalyst.incidencematgraph(trn.D_T)

    for lc in linkageclasses(trn)
        sg, _ = Graphs.induced_subgraph(g, lc)
        Graphs.is_strongly_connected(sg) || (trn.weaklyrev = false)
    end
    trn.checkedrev = true
    trn.weaklyrev
end

function linkageclasses(trn::Translation) 
    !isempty(trn.linkageclasses) && return trn.linkageclasses
    trn.linkageclasses = Graphs.connected_components(Catalyst.incidencematgraph(trn.D_T))
    trn.linkageclasses
end

function effectivedeficiency(trn::Translation) 
    trn.effectivedeficiency >= 0 && return trn.effectivedeficiency

    rn = trn.rn
    S = netstoichmat(rn)
    nc = size(trn.D_T, 1)
    l = length(linkageclasses(trn))
    trn.effectivedeficiency = nc - l - rank(S)
    trn.effectivedeficiency
end

function kineticdeficiency(trn::Translation) 
    trn.kineticdeficiency >= 0 && return trn.kineticdeficiency

    rn = trn.rn
    D = incidencemat(rn)
    nc = size(D, 1)
    l = length(linkageclasses(rn))

    trn.kineticdeficiency = nc - l - rank(trn.Y_T * D)
    trn.kineticdeficiency
end
