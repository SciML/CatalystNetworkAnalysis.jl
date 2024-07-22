const M::Float64 = 1E4
const ϵ::Float64 = 1E-1


# A network is discordant if there exists a nonzero σ ∈ image(S) and α ∈ ker(S) with the following sign properties: 
# 1) If α[r] != 0 for some reaction's index r, then the reaction's reactant complex must contain some species s for which sign(σ[s]) == sign(α[r])
# 2) If α[r] == 0 for some reaction r, either σ[s] == 0 for all s in the reactant complex, 
# or else there are two species s1, s2 in the reactant complex for which sign(σ[s1]) != sign(σ[s2])

function addconcordanceconstraints(model, rn::ReactionSystem) 
    α = model[:α]; σ = model[:σ]
    iszer = model[:σ_iszero]; ispos = model[:σ_ispos]; isneg = model[:σ_isneg]
    S = netstoichmat(rn); D = incidencemat(rn); Y = complexstoichmat(rn)
    s, r = size(S)

    @variable(model, allzero[1:r], Bin)
    @variable(model, allnonneg[1:r], Bin)
    @variable(model, allnonpos[1:r], Bin)
    @variable(model, bothposneg[1:r], Bin)
    @constraint(model, allzero + allnonneg + allnonpos + bothposneg == ones(r))

    for rxn in 1:r
        s = findfirst(==(-1), @view D[:, rxn])
        supp = findall(>(0), @view Y[:, s])
        numr = length(supp)

        # allzero ~ all(==(0), σ[supp]) <--> α[i] == 0        
        @constraint(model, sum(iszer[supp]) - M*(1 - allzero[rxn]) <= numr)
        @constraint(model, sum(iszer[supp]) + M*(1 - allzero[rxn]) >= numr)
        @constraint(model, α[rxn] - M*(1 - allzero[rxn]) <= 0)
        @constraint(model, α[rxn] + M*(1 - allzero[rxn]) >= 0)

        # allnonneg ~ all(>=(0), σ[supp]) <--> α[i] > 0
        @constraint(model, sum(isneg[supp]) - M*(1 - allnonneg[rxn]) <= 0)
        @constraint(model, sum(isneg[supp]) + M*(1 - allnonneg[rxn]) >= 0)
        @constraint(model, sum(ispos[supp]) + M*(1 - allnonneg[rxn]) >= 1)
        @constraint(model, α[rxn] + M*(1 - allnonneg[rxn]) >= ϵ)

        # allnonpos ~ all(<=(0), σ[supp]) <--> α[i] < 0
        @constraint(model, sum(ispos[supp]) - M*(1 - allnonpos[rxn]) <= 0)
        @constraint(model, sum(ispos[supp]) + M*(1 - allnonpos[rxn]) >= 0)
        @constraint(model, sum(isneg[supp]) + M*(1 - allnonpos[rxn]) >= 1)
        @constraint(model, α[rxn] - M*(1 - allnonpos[rxn]) <= -ϵ)

        # bothposneg ~ ∃s1, s2 σ[s1] > 0, σ[s2] < 0
        @constraint(model, sum(ispos[supp]) + M*(1 - bothposneg[rxn]) >= 1)
        @constraint(model, sum(isneg[supp]) + M*(1 - bothposneg[rxn]) >= 1)
    end
end

# Given a matrix S, return a linear programming model with the constraints that the
# desired vector μ will be sign-compatible with the image space of S. If the model
# already exists, add new sign compatibility constraints to the model.  

function signconstraintmodel(S::Matrix; model = nothing, var::String = "", in_subspace = false)
    (s, r) = size(S)

    model == nothing && begin 
        model = Model(HiGHS.Optimizer); 
        set_silent(model)
        set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)
        @objective(model, Min, 0)
    end
        
    coeffs = var*"_coeffs"
    model[Symbol(coeffs)] = @variable(model, [i = 1:r], base_name = coeffs) 
    model[Symbol(var)] = @variable(model, [i = 1:s], base_name = var)

    if in_subspace
        @constraint(model, S*model[Symbol(coeffs)] == model[Symbol(var)])
        return model
    end

    # Ensure that μ is sign-compatible with the stoichiometric subspace.
    ispos = var*"_ispos"; isneg = var*"_isneg"; iszer = var*"_iszero"

    model[Symbol(ispos)] = @variable(model, [i = 1:s], Bin, base_name = ispos)
    model[Symbol(isneg)] = @variable(model, [i = 1:s], Bin, base_name = isneg)
    model[Symbol(iszer)] = @variable(model, [i = 1:s], Bin, base_name = iszer) 

    @constraints(model, begin
       model[Symbol(iszer)] + model[Symbol(ispos)] + model[Symbol(isneg)] == ones(s) 
       sum(model[Symbol(iszer)]) <= s - 1 # Ensure that var is not the zero vector.

       # iszero = 1 --> var[i] == 0 <--> (S * coeffs)[i] == 0
       model[Symbol(var)] + M * (ones(s) - model[Symbol(iszer)]) ≥ zeros(s)
       model[Symbol(var)] - M * (ones(s) - model[Symbol(iszer)]) ≤ zeros(s)
       (S*model[Symbol(coeffs)]) + M * (ones(s) - model[Symbol(iszer)]) ≥ zeros(s)
       (S*model[Symbol(coeffs)]) - M * (ones(s) - model[Symbol(iszer)]) ≤ zeros(s)

       # isnegative = 1 --> var[i] < 0 <--> (S * coeffs)[i] < 0
       model[Symbol(var)] - M * (ones(s) - model[Symbol(isneg)]) ≤ -ones(s) * ϵ
       (S*model[Symbol(coeffs)]) - M * (ones(s) - model[Symbol(isneg)]) ≤ -ones(s) * ϵ

       # ispositive = 1 --> var[i] > 0 <--> (S * coeffs)[i] < 0
       model[Symbol(var)] + M * (ones(s) - model[Symbol(ispos)]) ≥ ones(s) * ϵ
       (S*model[Symbol(coeffs)]) + M * (ones(s) - model[Symbol(ispos)]) ≥ ones(s) * ϵ
    end)

    model
end


"""
    isconcordant(rn::ReactionSystem, atol=1e-12)

    Given a reaction network (and an absolute tolerance for the nullspace matrix below which entries should be zero), test whether the reaction network's graph has a property called concordance. A concordant network will not admit multiple equilibria in any stoichiometric compatibility class. The algorithm for this check follows Haixia Ji's PhD thesis, (Ji, 2011).  
"""

function isconcordant(rn::ReactionSystem) 
    S = netstoichmat(rn); (n, r) = size(S)

    numcols, kerS = nullspace_right_rational(ZZMatrix(S))
    kerS = Matrix{Int}(kerS[:, 1:numcols])

    model = signconstraintmodel(S, var = "σ") 
    signconstraintmodel(kerS, model = model, var = "α", in_subspace = true)
    addconcordanceconstraints(model, rn)

    optimize!(model)
    !is_solved_and_feasible(model)
end









###############################
###### OLD IMPLEMENTATION #####
###############################

# Given a sign pattern for σ, return the partial sign pattern for α and the indices that are free
function generate_α_signpattern(rn::ReactionSystem, σ::Vector) 
    S = netstoichmat(rn); (n, r) = size(S)

    # Free indices are ones that are unassigned, based on the sign pattern of σ. 
    α_sp = Int64[]; freeindices = Int64[]

    for rxn in 1:r 
        supp = findall(<(0), @view S[:, rxn])

        # To be a discordance, α must be 0 whenever all the species in the reactant complex are 0 in σ, 
        # and positive (negative) if they are all greater than (less than) or equal to 0. Otherwise, we impose no restrictions.  
        if all(==(0), σ[supp]) 
            push!(α_sp, 0)
        elseif all(>=(0), σ[supp]) 
            push!(α_sp, 1) 
        elseif all(<=(0), σ[supp]) 
            push!(α_sp, -1)
        else
            push!(freeindices, rxn)
        end
    end

    α_sp, freeindices
end

# Check if a vector v is sign-compatible with the image space of a matrix S
function issigncompatible(S::Matrix, v::Vector; freeindices::Vector{Int64}=Int[], NL = false)
    n, m = size(S)

    if length(v) + length(freeindices) != n
        error("The number of free signs and assigned signs does not sum to the length of the vector.") end

    assignedindices = deleteat!(collect(1:n), freeindices)
    model = Model(HiGHS.Optimizer)
    set_silent(model)

    # Determine which indices must be positive, negative, or zero in the solution to the LP. 
    zeroindices = assignedindices[findall(==(0), v)]
    negindices = assignedindices[findall(<(0), v)]
    posindices = assignedindices[findall(>(0), v)]

    @variable(model, coeffs[1:m])
    @objective(model, Min, 0)

    @constraint(model, (S*coeffs)[zeroindices] == zeros(length(zeroindices)))
    @constraint(model, (S*coeffs)[posindices] >= ones(length(posindices)))
    @constraint(model, (S*coeffs)[negindices] <= -ones(length(negindices)))

    optimize!(model)
    is_solved_and_feasible(model) ? true : false
end

function isconcordant_old(rn::ReactionSystem, atol=1e-12) 
    S = netstoichmat(rn); (n, r) = size(S)

    numcols, kerS = nullspace_right_rational(ZZMatrix(S))
    kerS = Matrix{Int}(kerS[:, 1:numcols])

    # Check whether there are species with fixed signs of zero. This is the case whenever the species does not appear in the support of any reaction vector. 
    fixedsigns = Int64[];
    for s in 1:n
        all(==(0), S[s, :]) && push!(fixedsigns, s)
    end
    sp = Int64[]

    # Initialize the model. 
    model = signconstraintmodel(S, var = "σ") 
    signconstraintmodel(kerS, model = model, var = "α", in_subspace = true)
    addconcordanceconstraints(model, S)
    
    # We check this by checking every possible sign pattern for σ by traversing a tree. 
    # Move forward until reaching a leaf node or an incompatible sign pattern for σ, then backtrack.  
   
    while true 
        println(sp)

        # Check whether the sign pattern for σ is compatible with image(S)
        if issigncompatible(S, sp, model, freeindices = collect(length(sp)+1:n))

            # If we have reached a leaf node of the tree, we have found a σ sign pattern compatible with image(S).
            if length(sp) == n

                # If we have reached the sign pattern of only zeros, we have checked every sign pattern. 
                # Since σ must be nonzero, there are no discordances and we return true for concordant. 
                all(==(0), sp) && return true

                # Each sign pattern for σ imposes restrictions on the allowable sign patterns for α, if (σ, α) is to be a discordance. 
                # We generate this α sign pattern here. 
                α_sp, freeidxs = generate_α_signpattern(rn, sp)

                # If the α sign pattern is compatible with ker(S), then we have found a discordance, and we return false for discordant. 
                issigncompatible(kerS, α_sp, freeindices = freeidxs) ? (return false) : (sp = movebackward(sp, fixedsigns))

            else
                # If we are not at a leaf node, and our partial sign pattern for σ is compatible with image(S), continue down the path. 
                sp = moveforward(sp, fixedsigns, n)
            end

        else
            # If our partial sign for σ is not compatible with image(S), backtrack and begin traversing a new branch. 
            sp = movebackward(sp, fixedsigns)
        end
    end
end

# Move forward simply adds + signs to the sign pattern, unless the species is a fixed sign, in which case we add 0. 

function moveforward(signpattern::Vector{Int64}, fixedsigns::Vector{Int64}, n::Int64) 
    while length(signpattern) + 1 ∈ fixedsigns
        push!(signpattern, 0)
    end

    length(signpattern) == n ? (return signpattern) : (return push!(signpattern, 1))
end

# Move backward triggers whenever we find a sign-pattern that is not compatible with image(S), pruning those branches of the tree

function movebackward(signpattern::Vector{Int64}, fixedsigns::Vector{Int64})
    if signpattern == [] return [] end

    # If we are currently at a fixed sign
    if length(signpattern) ∈ fixedsigns
        while length(signpattern) ∈ fixedsigns
            pop!(signpattern)
        end
        return movebackward(signpattern, fixedsigns)
    elseif signpattern[end] == 0
        pop!(signpattern)
        return movebackward(signpattern, fixedsigns)
    else
        # If we are currently at a positive or negative node, move to a different child of the mother node. 
        if signpattern[end] == 1
            pop!(signpattern)
            all(==(0), signpattern) ? push!(signpattern, 0) : push!(signpattern, -1)
            return signpattern
        else
            pop!(signpattern)
            push!(signpattern, 0)
            return signpattern
        end
    end
end
    

"""
    hasuniqueequilibria(rn::ReactionSystem)

    Check whether a reaction network has the capacity to admit multiple equilibria, for some choice of rate constants. 
"""

function hasuniqueequilibria(rn::ReactionSystem) 
    nps = get_networkproperties(rn)
    complexes, D = reactioncomplexes(rn)
    δ = deficiency(rn)
    concordant = isconcordant(rn)

    δ == 0 && return true 
    Catalyst.satisfiesdeficiencyone(rn) && return true 
    concordant && return true
    !concordant && ispositivelydependent(rn) && return false
    δ == 1 && return deficiencyonealgorithm(rn)
    
    error("The network is discordant and high deficiency, but this function currently cannot conclude whether the network has the potential to have multiple equilibria.")
end

# function isstronglyconcordant(rn::ReactionSystem) 
#     
# end
# 
# function isdegenerate(rn::ReactionSystem) 
#     
# end
# 
# function fullyopenextension(rn::ReactionSystem) 
#     
# end
# 
# function speciesreactiongraph(rn::ReactionSystem) 
#     s = numspecies(rn); sm = speciesmap(rn)
#     r = numreactions(rn); 
# 
#     G = Digraph(s+r)
#     adj = zeros(s+r, s+r)
#    
# end
