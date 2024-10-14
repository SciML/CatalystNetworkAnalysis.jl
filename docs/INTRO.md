### CatalystNetworkAnalysis.jl for Reaction Network Analysis
CatalystNetworkAnalysis.jl is a package that provides network analysis functionality
for reaction networks built using Catalyst.jl. The algorithms implemented by this 
package, along with their relevant citations, can be found in the [Algorithms] 
section of the documentation. 

Network analysis algorithms infer dynamical properties of reaction networks, such as 
the existence and uniqueness of steady states, from the graph structure of the network. 
This package implements algorithms that check for the existence of these graph structures. 

The algorithms currently implemented in the package draw from two broad approaches: 
deficiency theory (developed by Martin Feinberg et al.), which draws explicitly from the graph 
strcuture of the network, and approaches drawing from algebraic geometry, which analyzes the set 
of points in concentration space corresponding to the steady states, as well as the ideal that
corresponds to this set of points.
In these docs the algorithms used to determine reaction network properties will be cited - see the [Algorithms] page. 

The main function of the package is `networksummary(rn::ReactionSystem)`, which collects all the 
dynamical information that can be learned about the reaction network from its network structure. 

### Comparisons to existing packages
