using CatalystNetworkAnalysis
using SafeTestsets, Test
using UnPack, SBMLImporter, SBMLToolkit

@testset "CatalystNetworkAnalysis.jl" begin
    # Write your tests here.

    @time @safetestset "Concentration Robustness" begin include("ACR.jl")
    @time @safetestset "Concordance Helpers" begin include("concordancehelpers.jl") end
    @time @safetestset "Siphons" begin include("siphons.jl") end
    @time @safetestset "Persistence" begin include("persistence.jl") end
    @time @safetestset "Deficiency One Algorithm" begin include("deficiencyonealgorithm.jl") end
    @time @safetestset "Specific Stoichiometric Compatibility Class Functionality" begin include("specificscc.jl") end
end
