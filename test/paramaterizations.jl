include("./test_networks.jl")
let
    his_kinase = @reaction_network begin
        k1, X --> X_p
        (k2, k3), X_p + Y <--> X + Y_p
        k4, Y_p --> Y
    end
    true_parameterization = [
                            ]

    treeprods = matrixtree(his_kinase)

    @test false
end

let
    rn = EnvZ_OmpR

    true_parameterization = [
                            ]
    img = incidencematgraph(rn)
    distmx = Catalyst.adjacencymat(rn)
    treeprods(img, distmx)
    
    @test false
end

let
    rn = WNT

    true_parameterization = [
                            ]
    img = incidencematgraph(rn)
    distmx = Catalyst.adjacencymat(rn)
    treeprods(img, distmx)

    @test false
end
