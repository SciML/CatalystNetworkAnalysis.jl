"""
    ispersistent(rs::ReactionSystem)

    Checks if a reaction system is persistent, meaning that none of its species with positive concentration will go extinct (asymptotically approach 0). 
"""

function ispersistent(rs::ReactionSystem) 
    conslaws = conservationlaws(rs)
    conservative = !isempty(conservationlaws(rs))
    consistent = !isempty(cycles(rs))

    siphons = minimalsiphons(rs)
    conservative && consistent && all(!iscritical(s, conslaws), siphons)
end

"""
    minimalsiphons(rs::ReactionSystem)

    Constructs the set of minimal siphons of a reaction network, where a siphon is a set of species that can be "switched off," i.e. if the species each have concentration 0, the concentration of all the species will remain 0 for all time. A minimal siphon is one that does not contain a siphon as a subset.
"""

function minimalsiphons(rs::ReactionSystem, algorithm=:SMT) 
    if algorithm == :SMT 
        return minimalsiphons_smt(rs)
    else
        return minimalsiphons_alg(rs)
    end
end

# reimplement using Z3 - Satisfiability is a no-go (requires Z3 to be installed)

function minimalsiphons_smt(rs::ReactionSystem) 
    ns = numspecies(rs)
    sm = speciesmap(rs)
    @satvariable(specs[1:ns], Bool)
    constraints = [or(specs)] 
    siphons = Array{Int64}[]

    # Set up siphon constraints
    for rx in reactions(rs)
        subs = rx.substrates; prods = rx.products
        sub_idx = [sm[sub] for sub in subs]
        prod_idx = [sm[prod] for prod in prods]
        
        for p in prod_idx
            if isempty(subs)
                cons = not(specs[p])
            else
                cons = implies(specs[p], or([specs[s] for s in sub_idx]))
            end
            push!(constraints, cons)
        end
    end

    # Iterate until all siphons are found
    status = sat!(constraints..., solver=Z3())

    while status == :SAT
        siphon = findall(value(species))
        push!(siphons, siphon)
        push!(constraints, or(not.(species[siphon])))
        status = sat!(constraints..., solver=Z3())
    end

    return siphons
end

using Oscar

# Check if this can handle open reaction networks
function minimalsiphons_alg(rs::ReactionSystem) 
    sm = speciesmap(rs)
    specs = species(rs)
    rxns = reactions(rs)
    complexes, D = reactioncomplexes(rs)

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
        s = findfirst(==(-1), @view D[:,r])
        p = findfirst(==(1), @view D[:,r])
        polynomial = cm[s]*(cm[p] - cm[s])
        push!(ideal_generators, polynomial)
    end
    I = ideal(R, ideal_generators)

    siphons = [indexin(gens(prime), vars) for prime in minimal_primes(I)]
end

"""
    iscritical(s, conslaws)

    Checks if a siphon is critical, meaning that it does not contain the support of some conservation law. A reaction network with a critical siphon cannot be persistent.
"""
function iscritical(s::Vector, conslaws) 
    supports = [findall(!=(0), conslaws[i, :]) for i in 1:length(claws)]

    # If the support of any conservation law is contained in the siphon, then it is not critical
    all(sup->!issubset(sup, s), supports)
end
