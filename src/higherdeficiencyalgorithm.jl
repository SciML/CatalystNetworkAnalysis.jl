###################################
### HIGHER DEFICIENCY ALGORITHM ###
###################################

function higherdeficiencyalgorithm(rn::ReactionSystem) 
    # Choose an orientation. 
    S = netstoichmat(rn); rev_idxs = []
    for i in 1:size(S, 2)-1
        (S[:, i] == -S[:, i+1]) && push!(rev_idxs, i)
    end
    orientation_idxs = deleteat!(collect(1:size(S,2)), rev_idxs .+= 1)
    L_O = S[:, orientation_idxs];
    numcols, kerL_O = nullspace_right_rational(ZZMatrix(L_O)); 
    kerL_O = Matrix{Int64}(kerL_O)[:, 1:numcols]

    # Find the fundamental classes. 
    eq_classes = Dict()
    for i in 1:size(kerL_O, 1)
        row = kerL_O[i, :]; row = div.(row, gcd(row))
        r_i = orientation_idxs[i]
        haskey(eq_classes, row) ? push!(eq_classes[row], r_i) : eq_classes[row] = [r_i]
    end

    # Select representatives. 
    W_rev = Vector{Int}; W_irrev = Vector{Int}
    for rxs in values(eq_classes)
        irrev_rxs = filter!(r -> r ∉ rev_idxs, rxs)
        isempty(irrev_rxs) ? push!(W_rev, rxs[1]) : push!(W_irrev, irrev_rxs[1])
    end

    # Reorient. 

    # Build sign constraint model and add reaction constraints. 
    model = C.signconstraintmodel(S, model = model, var = "μ")
    C.signconstraintmodel(kerL_O, model = model, var = "g", in_subspace = true)
    C.signconstraintmodel(kerL_O, model = model, var = "h", in_subspace = true)
    
    prodmat = prodstoichmat(rn); submat = substoichmat(rn)
    addirrevconstraints(model, submat, W_irrev)
    addrevconstraints(model, submat, prodmat, W_rev)
end

# Irreversible Reactions: g > 0, h > 0, ρ = exp(y⋅μ)
function addirrevconstraints(model, prod, sub, idxs) 
    numrx = length(idx)
    @constraint(model, g[idxs] ≥ ϵ * ones(numrx))
    @constraint(model, h[idxs] ≥ ϵ * ones(numrx))

    expsubs = zeros(numrx)
    for i in idxs
        sub = @view submat[:, i]
        @constraint(model, g[i] == h[i]*exp(sub'*μ))
    end
end

# Reversible Reactions: 
# If g > 0, then one of the following holds:
#   ρ > exp(y⋅μ) > exp(y'⋅μ) 
#   ρ < exp(y⋅μ) < exp(y'⋅μ) 
#   ρ = exp(y⋅μ) = exp(y'⋅μ) 
#
# If g < 0, ρ > exp(y'⋅μ) > exp(y⋅μ)
#   ρ > exp(y'⋅μ) > exp(y⋅μ) 
#   ρ < exp(y'⋅μ) < exp(y⋅μ) 
#   ρ = exp(y'⋅μ) = exp(y⋅μ) 
#
# If g = 0, sgn(h) = sgn(exp(y⋅μ) - exp(y'⋅μ))
#
# Essentially we have nine possibilities. 

function addrevconstraint(model, prod, sub, idxs) 
    numrx = length(idxs)

    # sign of g: (pos, neg, zero)
    # relation of ρ and exp terms: (geq, eq, leq)
    varnames = ["posgeq", "poseq", "posleq", "neggeq", "negeq", "negleq", "zerogeq", "zeroeq", "zeroleq"]
    for var in varnames
        model[Symbol(var)] = @variable(model, [i=1:numrx], Bin, basename = var)
    end
    
    @constraints(model, begin
            sum([model[Symbol(var)] for var in varnames]) == ones(numrx)

            # posgeq == 1 --> g > 0 <--> ρ > exp(y⋅μ) > exp(y'⋅μ)
            model[Symbol("posgeq")] 

        end)
end

function findforestalbasis() 
    
end

function generateshelving() 
    
end
