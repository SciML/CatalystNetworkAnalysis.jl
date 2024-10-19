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
end
