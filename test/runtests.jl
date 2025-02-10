using CatalystNetworkAnalysis
using SafeTestsets, Test
using UnPack, SBMLImporter, SBMLToolkit

@testset "CatalystNetworkAnalysis.jl" begin
    # Write your tests here.

    @time @safetestset "Concentration Robustness" begin include("ACR.jl") end
    @time @safetestset "Concordance Helpers" begin include("concordancehelpers.jl") end
    @time @safetestset "Siphons" begin include("siphons.jl") end
    @time @safetestset "Persistence" begin include("persistence.jl") end
    @time @safetestset "Deficiency One Algorithm" begin include("deficiencyone.jl") end
    @time @safetestset "Specific Stoichiometric Compatibility Class Functionality" begin include("specificscc.jl") end
    @time @safetestset "Linear programming utilities" begin include("lp_utils.jl") end
    @time @safetestset "Network Translation" begin include("network_translation.jl") end

end
