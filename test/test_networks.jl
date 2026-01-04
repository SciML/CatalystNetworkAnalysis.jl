# Test networks (From Johnston et al, 2018)
using SBMLImporter
t = Catalyst.default_t()

MAPK = @reaction_network mapk begin
    (k1, k2), X + K <--> XK
    k3, XK --> Xp + K
    (k4, k5), Xp + K <--> XpK
    k6, XpK --> Xpp + K
    (k7, k8), Xpp + M <--> XppM
    k9, XppM --> XpM
    (k10, k11), XpM <--> Xp + M
    (k12, k13), Xp + M <--> Xp_M
    k14, Xp_M --> XM
    (k15, k16), XM --> X + M
end

zigzag = @reaction_network zigzag begin
    (k1, k2), X1 + X2 <--> X3
    k3, X3 --> X3 + X4
    k4, X4 --> ∅
    k5, X4 + X5 --> X5
    (k6, k7), X5 + X6 <--> X7
    k8, X8 --> X8 + X9
    k9, X9 --> X1 + X9
    k10, X9 --> X10 + X9
    k11, X1 --> ∅
    k12, X5 --> ∅
    k13, X9 --> ∅
    k14, X10 --> ∅
    k15, X7 + X9 --> X7
    k16, X4 + X9 --> X4
    (k17, k18), X10 + X11 <--> X12
    k19, X12 --> X5 + X11
    (k20, k21), X4 + X11 <--> X13
end

#####################################
### CONCENTRATION-ROBUST NETWORKS ###
#####################################

# From Johnston et al.
EnvZ_OmpR = @reaction_network EnvZ_OmpR begin
    (k1, k2), XD <--> X
    (k3, k4), X <--> XT
    k5, XT --> Xp
    (k6, k7), Xp + Y <--> XpY
    k8, XpY --> X + Yp
    (k9, k10), XT + Yp <--> XTYp
    k11, XTYp --> XT + Y
    (k12, k13), XD + Yp <--> XDYp
    k14, XDYp --> XD + Y
end

WNT = @reaction_network WNT begin
    (k1, k2), X1 <--> X2
    (k3, k4), X2 + X4 <--> X14
    k5, X14 --> X2 + X5
    (k6, k7), X5 + X8 <--> X16
    k8, X16 --> X4 + X8
    (k9, k10), X4 + X10 <--> X18
    k11, X18 --> X4
    (k12, k13), ∅ <--> X10
    (k14, k15), X3 + X6 <--> X15
    k16, X15 --> X3 + X7
    (k17, k18), X7 + X9 <--> X17
    k19, X17 --> X6 + X9
    (k20, k21), X6 + X11 <--> X19
    k22, X19 --> X6
    k23, X11 --> ∅
    (k24, k25), X11 + X12 <--> X13
    (k26, k27), X2 <--> X3
    (k28, k29), X5 <--> X7
    (k30, k31), X10 <--> X11
end

his_kinase = @reaction_network his_kinase begin
    k1, X --> Xp
    (k2, k3), Xp + Y <--> X + Yp
    k4, Yp --> Y
end

# From (S60) of Shinar et al., 2010. This network is deficiency two
feinberg_shinar_network = @reaction_network begin
    (k1, k2), XD <--> X
    (k3, k4), X <--> XT
    (k5), XT --> X_p
    (k6, k7), X_p + Y <--> X_pY
    k8, X_pY --> X + Y_p
    (k9, k10), XT + Y_p <--> XTY_p
    k11, XTY_p --> XT + Y
    (k12, k13), XD + Y_p <--> XDY_p
    k14, XDY_p --> XD + Y
end

acr_nets = [EnvZ_OmpR, WNT, his_kinase, feinberg_shinar_network]

###########################
### MULTIPLE EQUILIBRIA ###
###########################

# Some reaction networks with multiple equilibria, drawn from Gross et al, 2020
function oneSitePhosphorylation(n::Int64)
    rx = []
    n = 10
    @species S(t)[0:n] X(t)[1:n] Y(t)[1:n] E(t) F(t)
    @parameters k[1:n, 1:6]
    for i in 1:n
        push!(rx, Reaction(k[i, 1], [S[i - 1], E], [X[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i, 2], [X[i]], [S[i - 1], E], [1], [1, 1]))
        push!(rx, Reaction(k[i, 3], [S[i], F], [Y[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i, 4], [Y[i]], [S[i], F], [1], [1, 1]))
        push!(rx, Reaction(k[i, 5], [X[i]], [S[i], E], [1], [1, 1]))
        push!(rx, Reaction(k[i, 6], [Y[i]], [S[i - 1], F], [1], [1, 1]))
    end

    @named osp = ReactionSystem(rx, t, [S..., X..., Y..., E, F], vec(k))
    osp = complete(osp)
    return osp
end

function edelstein(n::Int64)
    rx = []
    n = 10
    @species A(t) B(t) D(t)[1:n]
    @parameters k[1:(4 * n + 2)]
    push!(rx, Reaction(k[1], [A], [A], [1], [2]))
    push!(rx, Reaction(k[2], [A], [A], [2], [1]))
    for i in 1:n
        i1, i2, i3, i4 = (4 * (i - 1) + 3, 4 * (i - 1) + 4, 4 * (i - 1) + 5, 4 * (i - 1) + 6)
        push!(rx, Reaction(k[i1], [A, B], [D[i]], [1, 1], [1]))
        push!(rx, Reaction(k[i2], [D[i]], [A, B], [1], [1, 1]))
        push!(rx, Reaction(k[i3], [D[i]], [B], [1], [1]))
        push!(rx, Reaction(k[i4], [B], [D[i]], [1], [1]))
    end
    @named edelstein = ReactionSystem(rx, t, [A, B, D...], collect(k))
    edelstein = complete(edelstein)
    return edelstein
end

function cellDeathNetwork(n::Int64)
    rx = []
    n = 10
    @species X(t) Y(t)
    @parameters k[1:n, 1:n]
    for i in 1:n
        for j in 1:n
            i >= j && continue
            push!(rx, Reaction(k[i, j], [X, Y], [X, Y], [n - i, i], [n - j, j]))
        end
    end
    @named rs = ReactionSystem(rx, t)
    rs = complete(rs)
    return rs
end

#################################
### MASS-ACTION SBML NETWORKS ###
#################################

netdir = joinpath(@__DIR__, "mass-action-networks/")
ma_nets = Dict()

for file in readdir(netdir)
    name, ext = splitext(file)
    path = joinpath(netdir, file)
    try
        prn, cb = load_SBML(path)
        ma_nets[name] = prn
    catch e
        continue
    end
end
