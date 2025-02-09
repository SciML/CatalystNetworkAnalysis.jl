const C = CatalystNetworkAnalysis 

"""
    symbolic_steady_states(rn::ReactionSystem)

    Given a reaction system, compute the parameterization of the system's steady state in terms of the parameters.
    1) For generalized mass action systems with kinetic deficiency zero. 
    2) Algorithm for positive kinetic deficiency systems. 
    3) Monomial parameterization for systems with toric steady states (Millan et al. 2012)

    The output of this function is a symbolic function representing the parameterization. It takes some subset of the species and maps them into a steady state.  

    This function will return `nothing` if no algorithm currently is implemented to find symbolic steady states for the case.
"""
function symbolic_steady_states(rn::ReactionSystem; variable_search = false) 
    trn = WRDZ_translation(rn)

    if !isnothing(trn) && kineticdeficiency(trn) == 0
        return sss_wrdz(trn)
    elseif !isnothing(trn) && kineticdeficiency(rn) > 0
        return sss_pos_δ(rn)
    elseif istoric(rn)
        return sss_toric(rn)
    else
        return nothing
    end

    # κ ∘ τ
    # Find a spanning forest
    # Find B such that im B = ker M^T
    # H is a generalized inverse of M
end

# Monomial parameterizations of systems that are WRDZ
function sss_wrdz(rn::ReactionSystem) 
    translation = WRDZ_translation(rn)

    img = incidencematgraph(translation.rn)
    D = incidencemat(translation.rn)

    # Get the spanning forest
    forest = vcat ∘ map(Graphs.connected_components(img)) do cc 
        sg, vmap = induced_subgraph(img, cc)
        tree = kruskal_mst(sg)
        edges = map(e -> Edge(vmap[src(e)], vmap[dst(e)]), tree)
    end
    edgeidxs = indexin(forest, collect(edges(img)))

    M = translation.Y_T * (@view D[:, edgeidxs])
    distmx = Catalyst.ratematrix(translation.rn)
    treeprods = matrixtree(img, distmx)

    n, ϵ = size(M)
    s = rank(Y_T*D)

    ncol, B = nullspace_right_rational(ZZMatrix(M'))
    B = @view B[:, 1:ncol]
    κ = map(forest) do e
        treeprods[src(e)] / treeprods[dst(e)]
    end
    @variables τ[1:n-s]

    H = pinv(M')
    C.matrixpower(κ, H') .* C.matrixpower(τ, B')
end

function sss_pos_δ(rn::ReactionSystem) 
    
end

function toric(rn::ReactionSystem) 
    treeproducts = matrixtree(rn)
    spanforest
end
