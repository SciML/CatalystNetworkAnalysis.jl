# Sign conditions for injectivity, following
#
# (1) Müller, S.; Feliu, E.; Regensburger, G.; Conradi, C.; Shiu, A.; Dickenstein, A. Sign Conditions for Injectivity of Generalized Polynomial Maps with Applications to Chemical Reaction Networks and Real Algebraic Geometry. arXiv October 30, 2014. http://arxiv.org/abs/1311.5493 (accessed 2024-08-15).

# is there a vector x such that: 
#   x is sign-compatible with ker(S)
#   y is sign-compatible with S^*
#   x is sign-compatible with V(y)

"""
    isinjective(rn::ReactionSystem)

    Determine whether a reaction network system is injective, i.e. will have at most one steady state.
"""
function isinjective(rn::ReactionSystem)
    S = netstoichmat(rn)
    Y = complexstoichmat(rn)
    D = incidencemat(rn)
    kerS = nullspace_right_rational(S)

    # Kinetic order matrix
    r, s = (length(reactions(rn)), length(species(rn)))
    V = zeros(Int64, r, s)

    for i in 1:length(reactions(rn))
        src = findfirst(==(-1), @view D[:, i])
        V[i, :] .= Y[src, :]
    end

    model = add_sign_constraints(kerS, var_name = "x")
    add_sign_constraints(S; model, var_name = "y")

    # Add injectivity constraints.
    x = model[:x]
    y = model[:y]
    iszer = model[:x_iszero]
    ispos = model[:x_ispos]
    isneg = model[:x_isneg]

    # Ensure that V*y, x have the same sign pattern
    #   V*y > 0 <--> x > 0
    #   V*y == 0 <--> x == 0
    #   V*y < 0 <--> x < 0
    @constraints(model, begin
        # iszero = 1 --> x[i] == 0 <--> (V * y)[i] == 0
        x + M * (ones(s) - iszer) ≥ zeros(s)
        x - M * (ones(s) - iszer) ≤ zeros(s)
        V*y + M * (ones(s) - iszer) ≥ zeros(s)
        V*y - M * (ones(s) - iszer) ≤ zeros(s)

        # isnegative = 1 --> x[i] < 0 <--> (V * y)[i] < 0
        x - M * (ones(s) - isneg) ≤ -ones(s) * ϵ
        V*y - M * (ones(s) - isneg) ≤ -ones(s) * ϵ

        # ispositive = 1 --> x[i] > 0 <--> (V * y)[i] > 0
        x + M * (ones(s) - ispos) ≥ ones(s) * ϵ
        V*y + M * (ones(s) - ispos) ≥ ones(s) * ϵ
    end)

    optimize!(model)
    !is_solved_and_feasible(model)
end
