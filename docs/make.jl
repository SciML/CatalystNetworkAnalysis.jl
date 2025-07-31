using Documenter
using CatalystNetworkAnalysis, Catalyst

makedocs(
    sitename = "CatalystNetworkAnalysis.jl",
    authors = "Vincent Du",
    format = Documenter.HTML(prettyurls = (get(ENV, "CI", nothing) == "true")),
    modules = [CatalystNetworkAnalysis, Catalyst],
    doctest = false,
    clean = true,
    pages = Any[
        "Home" => "index.md",
        "Network Analysis Algorithms" => "Algorithms.md",
        "Roadmap" => "ROADMAP.md"
    ],
    warnonly = [:missing_docs]
)

deploydocs(
    repo = "github.com/SciML/CatalystNetworkAnalysis.jl.git";
    push_preview = true
)
