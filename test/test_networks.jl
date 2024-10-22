# Test networks (From Johnston et al, 2018)

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


# From Johnston et al 
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
