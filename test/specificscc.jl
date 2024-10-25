t = Catalyst.default_t()

# Test the mixed volume function
let
    rn = @reaction_network begin
        k1, A --> 0
        k2, A --> 2A
        (k3, k4), 2A <--> 3A
    end
    u0 = []

    modified_SFR(rn, u0)

    # Programmatically generating a cell death reaction network
    t = Catalyst.default_t()
    rx = []; n = 10
    @species X(t) Y(t)
    @parameters k[1:n, 1:n] = zeros(n, n)
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

    # Edelstein network

end
