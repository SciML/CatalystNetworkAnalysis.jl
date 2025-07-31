### CatalystNetworkAnalysis.jl for Reaction Network Analysis

CatalystNetworkAnalysis.jl is a package that provides network analysis functionality for reaction networks built using Catalyst.jl. The algorithms implemented by this package, along with their relevant citations, can be found in the [Network Analysis Algorithms Overview](@ref)
section of the documentation.

Network analysis algorithms infer dynamical properties of reaction networks, such as the existence and uniqueness of steady states, from the graph structure of the network and the structure of their ODEs. This package implements algorithms that take input reaction networks, and perform this analysis. The benefit of this approach is that key properties of the network may be inferred *without* simulating them, and we can infer properties that are true of the network regardless of its specific choice of reaction rate constants or initial conditions.

The algorithms currently implemented in the package draw from two broad approaches: deficiency theory (developed by Martin Feinberg et al.), which draws explicitly from the graph structure of the network, and approaches drawing from algebraic geometry, which analyzes the set of points in concentration space corresponding to the steady states, as well as the ideal that corresponds to this set of points. More information about the properties, and the algorithms implemented in this package  can be found in the [Algorithms] page.

### Network Summary

The main function of the package is `networksummary(rn::ReactionSystem)`, which collects all the dynamical information that can be learned about the reaction network from its network structure into a readable format. To see how it works, let us consider a simple reaction network that displays many curious dynamical properties called the Edelstein network:

```@example intro
using CatalystNetworkAnalysis, Catalyst
edelstein = @reaction_network begin
    (k1, k2), A <--> 2A
    (k3, k4), A + B <--> C
    (k5, k6), C <--> B
end
```

The network summary for this reaction network shows the following:

```@example intro
networksummary(edelstein)
```

### Specifying Parameter Values and Initial Conditions

Many of the algorithms implemented in this package will optionally take parameter assignments `p` or initial conditions `u0`. In some cases, these will be required (e.g. complex balance will require parameters, and mixed volume will require initial conditions). If `p` or `u0` are optional, they will be keyword functions in arguments, and if they are required for a function, they will be a positional argument. The `networksummary` function will likely have different outputs depending on whether `u0` or `p` or both are provided.

These will only be accepted in the standardized format of a `Dict`, a Vector of `Pair`, or a Tuple of `Pair`. The keys of the pairs should be either symbols or Symbolics. For example, to provide parameters for the Edelstein network above, any of the following would work:

```@example intro
pmap = [:k1 => 1.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.0, :k5 => 1.0, :k6 => 1.0]
pmap = (:k1 => 1.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.0, :k5 => 1.0, :k6 => 1.0)
pmap = Dict([:k1 => 1.0, :k2 => 1.0, :k3 => 1.0, :k4 => 1.0, :k5 => 1.0, :k6 => 1.0])
pmap = Dict(zip(parameters(edelstein), ones(6))) # Creates dict from symbolic pairs
```

In order to convert from a symbol dictionary to a Symbolics dictionary, one can use the `symmap_to_varmap` function from Catalyst.

### Comparisons to existing packages

A variety of packages for network analysis exist, including the following:

 1. [Chemical Reaction Network Toolbox](https://zenodo.org/records/5149266) (Feinberg et al.)
 2. [CoNTRoL](https://control.math.wvu.edu/) (Johnston et al.)

### Starting with Julia and Catalyst

Please see the [starting with Julia](https://docs.sciml.ai/Catalyst/stable/introduction_to_catalyst/catalyst_for_new_julia_users/) section of Catalyst.jl's documentation for tips on getting started with Julia in general.

If you have a reaction network stored in some other file format, such as `.net`, SBML, and others, please reference the [loading chemical reaction networks](https://docs.sciml.ai/Catalyst/stable/model_creation/model_file_loading_and_export/) tutorial in the Catalyst documentation. There are sections for loading `.net` and SBML files. Once the reaction network is loaded as a Catalyst `ReactionSystem`, it can be passed to the network analysis algorithms defined in this function.
