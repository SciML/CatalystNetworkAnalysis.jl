import CatalystNetworkAnalysis as C
using Catalyst
using Graphs

##########################################################
### Test weakly reversible + deficiency zero algorithm ###
##########################################################

# Kinase example
let
    kinase = @reaction_network begin
        r1, X --> Xp
        (r2, r3), Xp + Y <--> X + Yp
        r4, Yp --> Y
    end

    translation = C.WRDZ_translation(kinase)
    @test translation.Y_T * translation.D_T == netstoichmat(kinase)
    @test C.deficiency(translation) == 0
    @test C.isweaklyreversible(translation)
end

let
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

    # Test common source reactions are found correctly
    S = netstoichmat(zigzag)
    csr = C.common_source_reactions(zigzag)
    csr_true = [[2,3],
                [9,10,13],
                [18,19]]
    @test issetequal(csr, csr_true)

    # Test elementary flux modes 
    efms = C.elementary_flux_modes(zigzag)
    @test all(col -> iszero(S*col), eachcol(efms))

    efm_supports = [findall(>(0), efm) for efm in eachcol(efms)]
    efm_parts = C.efm_partitions(zigzag, efm_supports, csr)
    efm_true = [[1,2,3,4,5],
                [6,7],
                [8,9,10,11,12,13,14,15,16,17,18,19],
                [20,21]]

    @test all(issetequal.(efm_parts, efm_true))

    rrg_mat = C.construct_rr_graph(zigzag)

    # CS constraints
    for rs in csr 
        edges = rrg_mat[:, rs[1]]
        @test all(==(edges), eachcol(rrg_mat[:, rs]))
    end

    # Partition constraints
    ccs = connected_components(SimpleDiGraph(rrg_mat))
    @test all(issetequal.(ccs, efm_parts)) 
    
    # EFM constraints
    for supp in efm_supports
        l = length(supp)
        part = @view rrg_mat[supp, supp]
        @test sum(part, dims=2) == ones(Bool, l, 1)
        @test sum(part, dims=1) == ones(Bool, 1, l)
        @test length(connected_components(Graphs.SimpleDiGraph(part))) == 1
    end

    translation = C.WRDZ_translation(kinase)
    @test translation.Y_T * translation.D_T == netstoichmat(kinase)
    @test C.deficiency(translation) == 0
    @test C.isweaklyreversible(translation)
end

let
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
        (k15, k16), XM <--> X + M
    end

    S = netstoichmat(MAPK)
    csr = C.common_source_reactions(MAPK)
    csr_true = [[2,3],
                [5,6],
                [8,9],
                [11,12],
                [13,14]]
    @test all(issetequal.(csr, csr_true))

    efms = C.elementary_flux_modes(MAPK)
    @test all(col -> iszero(S*col), eachcol(efms))

    efm_supports = [findall(>(0), efm) for efm in eachcol(efms)]
    efm_parts = C.efm_partitions(MAPK, efm_supports, csr)
    @test issetequal(efm_parts[1], collect(1:length(reactions(MAPK))))

    rrg_mat = C.construct_rr_graph(MAPK)

    # CS constraints
    for rs in csr 
        edges = rrg_mat[:, rs[1]]
        @test all(==(edges), eachcol(rrg_mat[:, rs]))
    end

    # Partition constraints
    ccs = connected_components(SimpleDiGraph(rrg_mat))
    @test all(issetequal.(ccs, efm_parts)) 
    
    # EFM constraints
    for supp in efm_supports
        l = length(supp)
        part = @view rrg_mat[supp, supp]
        @test sum(part, dims=2) == ones(Bool, l, 1)
        @test sum(part, dims=1) == ones(Bool, 1, l)
        @test length(connected_components(Graphs.SimpleDiGraph(part))) == 1
    end

    translation = C.WRDZ_translation(kinase)
    @test translation.Y_T * translation.D_T == netstoichmat(kinase)
    @test C.deficiency(translation) == 0
    @test C.isweaklyreversible(translation)
end

## Test generation of steady-state parameterizations
let
    rn = @reaction_network begin
        (k1, k2), A + B <--> C
        k3, C --> 2A
        k4, 2A --> A + B
        (k5, k6), A <--> B
    end

    kinase = @reaction_network begin
        r1, X --> Xp
        (r2, r3), Xp + Y <--> X + Yp
        r4, Yp --> Y
    end
end
