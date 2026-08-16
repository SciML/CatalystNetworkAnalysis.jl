### TODO
#   endotactic networks
#   permanence in networks
#   improve efficiency of siphon detection

"""
    ispersistent(rs::ReactionSystem)

Determine whether a reaction system is persistent.

A persistent system does not allow a species that starts with positive
concentration to asymptotically approach zero. The current structural test is
conclusive for some conservative systems and returns `:INCONCLUSIVE` for the
remaining cases.

# Arguments

- `rs`: reaction system to analyze.

# Returns

One of `:PERSISTENT`, `:NOT_PERSISTENT`, or `:INCONCLUSIVE`.
"""
function ispersistent(rs::ReactionSystem)
    siphons = minimalsiphons(rs)
    conservative = isconservative(rs)
    consistent = isconsistent(rs)
    S = netstoichmat(rs)

    # Conservative case
    if conservative
        all(s -> !iscritical(s, S), siphons) && return :PERSISTENT
        !consistent && return :NOT_PERSISTENT
    end

    return :INCONCLUSIVE
end

###############
### SIPHONS ###
###############

"""
    minimalsiphons(rs::ReactionSystem)

Construct the minimal siphons of a reaction system.

A siphon is a set of species that can be switched off: if every species in
the set has zero concentration, the species in the set remain zero. A minimal
siphon contains no siphon as a strict subset.

# Arguments

- `rs`: reaction system to analyze.

# Keyword Arguments

- `algorithm`: either `:SMT` (the default) for the satisfiability-based
  method or `:ALG` for the algebraic method.

# Returns

A vector of vectors of one-based species indices, one vector for each minimal
siphon.

# Throws

An error if `algorithm` is neither `:SMT` nor `:ALG`.
"""
function minimalsiphons(rs::ReactionSystem; algorithm = :SMT)
    if algorithm == :SMT
        return minimalsiphons_smt(rs)
    elseif algorithm == :ALG
        return minimalsiphons_alg(rs)
    else
        error("Invalid algorithm specified")
    end
end

function minimalsiphons_smt(rs::ReactionSystem)
    ns = numspecies(rs)
    sm = speciesmap(rs)

    # We encode the problem as a Boolean satisfiability problem. In a siphon search, species that belong to the siphon have a value of 1, and those that do not have a value of 0
    @satvariable(specs[1:ns], Bool)

    # Our initial constraint requires that there is at least one element in the siphon.
    constraints = [or(specs)]
    siphons = Array{Int}[]

    # Each reaction adds some constraint to our satisfiability problem.
    for rx in reactions(rs)
        # Determine substrate and product species for the given reaction.
        subs = rx.substrates
        prods = rx.products
        sub_idx = [sm[sub] for sub in subs]
        prod_idx = [sm[prod] for prod in prods]

        # Add constraints as such:
        # If the reaction has ∅ as a substrate complex, then it cannot be a member of a siphon.
        # If s is produced by the reaction, then s = 1 implies that there is some species in the substrate complex that is also equal to 1.

        for p in prod_idx
            if isempty(subs)
                cons = not(specs[p])
            else
                cons = implies(specs[p], or([specs[s] for s in sub_idx]))
            end
            push!(constraints, cons)
        end
    end

    # Solve the CSP to find a siphon.
    status = sat!(constraints..., solver = Z3())

    # Any time we find a siphon, we must add another constraint in order to ensure that the siphons are minimal. To disallow
    while status == :SAT
        siphon = findall(Satisfiability.value(specs))
        push!(siphons, siphon)
        push!(constraints, or(not.(specs[siphon])))
        status = sat!(constraints..., solver = Z3())
    end

    return removesupersets(siphons)
end

function removesupersets(indexsets)
    indexsets = sort(indexsets, by = length)
    minimalsets = Array{Int64}[]

    for s in indexsets
        if !any(ms -> issubset(ms, s), minimalsets)
            push!(minimalsets, s)
        end
    end
    return minimalsets
end

# TODO: Check if this can handle open reaction networks
function minimalsiphons_alg(rs::ReactionSystem)
    sm = speciesmap(rs)
    specs = species(rs)
    complexes, D = reactioncomplexes(rs)
    rxns = reactions(rs)

    R, vars = polynomial_ring(QQ, string.(species(rs)))

    cm = []
    for c in complexes
        if isempty(c)
            push!(cm, 0)
        else
            monomial = prod([vars[rce.speciesid]^rce.speciesstoich for rce in c])
            push!(cm, monomial)
        end
    end

    ideal_generators = []
    for r in 1:length(rxns)
        s = findfirst(==(-1), @view D[:, r])
        p = findfirst(==(1), @view D[:, r])
        polynomial = cm[s] * (cm[p] - cm[s])
        push!(ideal_generators, polynomial)
    end
    I = ideal(R, ideal_generators)

    return siphons = [indexin(gens(prime), vars) for prime in minimal_primes(I)]
end

"""
    iscritical(siphon, S)

Check whether a siphon is critical.

A siphon is critical when it does not contain the support of a positive
conservation law. A reaction network with a critical siphon cannot be
persistent.

# Arguments

- `siphon`: one-based species indices defining the siphon.
- `S`: stoichiometric matrix for the reaction system.

# Returns

`true` when `siphon` is critical and `false` otherwise.
"""
function iscritical(siphon::Vector, S::Matrix)
    # Takes the rows of the stoichiometric matrix corresponding to the siphon species
    S_r = S[siphon, :]

    # If there is a non-negative vector in the nullspace of S_red', then there is a positive conservation law with a support that is the subset of the siphon, and the siphon is not critical
    conslaws_r = conservationlaws(S_r)
    return !has_positive_solution(copy(conslaws_r'), nonneg = true)
end
