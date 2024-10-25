# CatalystNetworkAnalysis

[![Build Status](https://github.com/vyudu/CatalystNetworkAnalysis.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/vyudu/CatalystNetworkAnalysis.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![ColPrac: Contributor's Guide on Collaborative Practices for Community Packages](https://img.shields.io/badge/ColPrac-Contributor's%20Guide-blueviolet)](https://github.com/SciML/ColPrac)
[![SciML Code Style](https://img.shields.io/static/v1?label=code%20style&message=SciML&color=9558b2&labelColor=389826)](https://github.com/SciML/SciMLStyle)

CatalystNetworkAnalysis is a package providing analysis utilities for reaction networks defined in [Catalyst.jl](https://docs.sciml.ai/Catalyst/stable/). CatalystNetworkAnalysis provides the ability to determine dynamic properties of reaction networks, such as the ability to admit multiple steady states, the existence of species with concentrations robust to perturbation, persistence, and more. The algorithms implemented in the package extract information from the graph structure of `ReactionSystems` or the algebraic structure of their polynomial differential equations. For a full set of dynamical properties of networks that this package can shed light on, along with citations for the algorithms implemented, please see the [documentation].  

### Examples
To use this package, first install it and Catalyst from the Julia package manager: 
```julia
pkg> add Catalyst, CatalystNetworkAnalysis
```

Then, import a desired reaction network for analysis using tools like [ReactionNetworkImporters](https://github.com/SciML/ReactionNetworkImporters.jl) or create a chemical reaction network, e.g. by using Catalyst's DSL. Note that the network must be of Catalyst's `ReactionSystem`. 
```julia
# The Horn-Jackson network
rn = @reaction_network begin
    k1, 3A --> 2A + B
    k2, 2A + B --> 3B
    k1, 3B --> 2A + B
    k2, 2A + B --> 3A
end
```

Finally we can pass the reaction network into CatalystNetworkAnalysis's utilities. The main function is a `networksummary` that collates dynamical information computed for the network at-a-glance.
```julia
networksummary(rn)
```

Some functions will optionally take initial conditions `u0` or parameter maps `p`, which should be specified as `Dict`, `Vector{Pair}`, or `Tuple{Pair}`. 
```julia
p = [:k1 => 1., k2 => 0.15]
networksummary(rn; p = p)
```

### Getting help or getting involved

Catalyst and CatalystNetworkAnalysis developers are active on the [Julia Discourse](https://discourse.julialang.org/) and the [Julia Slack](https://julialang.slack.com/) channels #sciml-bridged and #sciml-sysbio. For bugs or feature requests, [open an issue](https://github.com/SciML/CatalystNetworkAnalysis.jl/issues).
