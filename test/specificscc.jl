t = Catalyst.default_t()

# Test the mixed volume function. These test networks are drawn from Gross, Hill 2020. 
let
    rn = @reaction_network begin
        k1, A --> 0
        k2, A --> 2A
        (k3, k4), 2A <--> 3A
    end

    modified_SFR(rn, u0)

    # Programmatically generating a cell death reaction network
    t = Catalyst.default_t()
    rx = []; n = 10
    @species X(t) Y(t)
    @parameters k[1:n, 1:n] 
    for i in 1:n
        for j in 1:n
            i >= j && continue
            push!(rx, Reaction(k[i, j], [X, Y], [X, Y], [n-i, i], [n-j, j]))
        end
    end
    @named rs = ReactionSystem(rx, t)
    rs = complete(rs)

    u0 = [:X => 1., :Y => 2.]
    @test C.mixedvolume(rs, u0) == n - 2
end

# Edelstein network
let
    rx = []; n = 10
    @species A(t) B(t) D(t)[1:n]
    @parameters k[1:(4*n+2)]
    push!(rx, Reaction(k[1], [A], [A], [1], [2]))
    push!(rx, Reaction(k[2], [A], [A], [2], [1]))
    for i in 1:n
        i1, i2, i3, i4 = (4*(i-1)+3, 4*(i-1)+4, 4*(i-1)+5, 4*(i-1)+6)
        push!(rx, Reaction(k[i1], [A, B], [D[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i2], [D[i]], [A, B], [1], [1, 1]))
        push!(rx, Reaction(k[i3], [D[i]], [B], [1], [1]))
        push!(rx, Reaction(k[i4], [B], [D[i]], [1], [1]))
    end
    @named edelstein = ReactionSystem(rx, t, [A, B, D...], collect(k))
    edelstein = complete(edelstein)
    u0 = Dict(zip(species(edelstein), ones(n+2)))
    @test C.mixedvolume(edelstein, u0) == 3
end

# One-site phosphorylation network
let
    rx = []; n = 10
    @species S(t)[0:n] X(t)[1:n] Y(t)[1:n] E(t) F(t)
    @parameters k[1:n, 1:6] 
    for i in 1:n
        push!(rx, Reaction(k[i, 1], [S[i-1], E], [X[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i, 2], [X[i]], [S[i-1], E], [1], [1, 1]))
        push!(rx, Reaction(k[i, 3], [S[i], F], [Y[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i, 4], [Y[i]], [S[i], F], [1], [1, 1]))
        push!(rx, Reaction(k[i, 5], [X[i]], [S[i], E], [1], [1, 1]))
        push!(rx, Reaction(k[i, 6], [Y[i]], [S[i-1], F], [1], [1, 1]))
    end

    @named osp = ReactionSystem(rx, t, [S..., X..., Y..., E, F], vec(k))
    osp = complete(osp)
    u0 = Dict(zip(species(osp), ones(3*n+3)))
    # @test C.mixedvolume(osp, u0) == div((n+1)*(n+4), 2) - 1
end
