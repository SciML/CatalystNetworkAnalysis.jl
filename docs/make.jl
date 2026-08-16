using Documenter
using CatalystNetworkAnalysis, Catalyst

makedocs(
    sitename = "CatalystNetworkAnalysis.jl",
    authors = "Vincent Du",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [CatalystNetworkAnalysis],
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "Network Analysis Algorithms" => "Algorithms.md",
        "Roadmap" => "ROADMAP.md",
    ]
)

deploydocs(
    repo = "github.com/SciML/CatalystNetworkAnalysis.jl.git";
    push_preview = true
)
