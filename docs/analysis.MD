# The Goal of Network Analysis

Network analysis algorithms infer dynamical properties of reaction networks, such as 
the existence and uniqueness of steady states, from the graph structure of the network. 
This package implements algorithms that check for the existence of these graph structures. 

The algorithms currently implemented in the package draw from two broad approaches: 
deficiency theory (developed by Martin Feinberg et al.), which draws explicitly from the graph 
strcuture of the network, and approaches drawing from algebraic geometry, which analyzes the set 
of points in concentration space corresponding to the steady states, as well as the ideal that
corresponds to this set of points.
In these docs the algorithms used to determine reaction network properties will be cited. 

The main function of the package is `networksummary(rn::ReactionSystem)`, which collects all the 
dynamical information that can be learned about the reaction network from its network structure. 

There are a few dynamical properties of particular interest that network approaches can shed
some light on. 
1. Existence of steady states
2. Uniqueness of steady states
3. Persistence
4. Concentration Robustness

The dynamics of a chemical reaction network depend on the forms of the rate laws and the parameters
chosen, particularly the rate constants. A number of the results in this package are dependent
on the choice of rate constants. We can distinguish between a chemical reaction network and a
kinetic system, which is a chemical reaction network along with an assignment of parameters. 

The following checks work for chemical reaction networks broadly: 
1. .

And the following checks are meant for kinetic systems (i.e. a set of parameters must be supplied): 
1. 

A final key notion is that of the **stoichiometric compatibility class**, which is the set of points 
in concentration space that are reachable from a given initial condition. When we say that a network
is multistationary, or has the capacity for multiple steady states, what that means is that there is
*some* stoichiometric compatibility class for which this is true. This means that there may be only one
or zero steady states in some other stoichiometric compatibility class. 

### Existence of Steady States


### Uniqueness of Steady States
This functionality is in the function `hasuniqueequilibria(rn)`.

The uniqueness of steady states checks proceeds as follows. The function terminates every time a
network is conclusively determined to have unique steady states or multi-stationarity.  

1. **Deficiency one and deficiency zero theorems**: a reaction network that satisfies the conditions of 
these two theorems will be guaranteed to have a unique steady state for every stoichiometric
compatibility class. 
2. **Deficiency one and higher deficiency algorithm**: Both of these algorithms construct inequality systems
that, if satsifiable, would indicate that a reaction network has the ability to admit multiple steady
states in a given stoichiometric compatibility class. 
3. **Concordance**: This is a graph property that guarantees the uniqueness of steady states, since it
guarantees that the species-formation rate function is injective. As a result, there is only one solution
to $\frac{dc}{dt} = 0$. 
4. **Monotone systems theory**: 
5. **Degree theory**: 
6. **Embedded multi-stationary networks**: 

Additionally, one might check the functions `iscomplexbalanced()` and `isdetailedbalanced()` implemented
in the main Catalyst package, which work for kinetic systems. Both of these conditions are sufficient 
to guarantee that the reaction network will have unique steady states, *for that choice of parameters*,
and that these steady states will be *asymptotically stable* within each stoichiometric compatibility class.

Eventually the function will have further tooling to determine properties of the steady states. 

### Number of Steady States
Upper Bound
It is known that the number of steady states that a polynomial ODE system can admit is bounded 
above by a quantity called the [mixed volume](). 

Bistability
Bistability is a very important property of reaction networks, since switchlike behavior in biological 
networks is realized through bistability. The network analysis package. 

For further stability analysis, please see [this section]() of the main Catalyst documentation.

### Persistence
A chemical reaction network is persistent if each species does not "die out" - that is, if the
initial condition of the network has positive concentration in a given species, that species'
concentration will never go to zero. 

It is not currently possible to conclusively determine persistence in every case. The approach
used by this package checks the necessary and sufficient conditions outlined in this paper. Any
persistent network will be. Note that it is currently a conjecture that any weakly reversible
network will be persistent. 

### Concentration Robustness
A species in a reaction network has concentration robustness if its concentration in a given
steady state of the reaction network is the same. 

1. Deficiency results: these are very limited to 
2. Algebraic results: 
