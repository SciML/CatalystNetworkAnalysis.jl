# Network Analysis Algorithms Overview

There are a few dynamical properties of particular interest that network approaches can shed some light on.

 1. Existence + uniqueness of steady states
 2. Persistence
 3. Concentration Robustness

Algorithms used by this package to decide these properties are detailed on this page. A significant number of the results can be found in the *Foundations of Chemical Reaction Network Theory* textbook, by Martin Feinberg[^1].

In general these properties will depend on the choices of parameters for the system and on the initial conditions (which can be provided as a dictionary, vector, or tuple of `Pair`, as described in the homepage). Many networks have multiple steady states for some set of rate constants but not others.

Additionally, the dynamical properties of the system may depend on the initial conditions of the network, which define the **stoichiometric compatibility class** that the network evolves in. The stoichiometric compatibility class is the set of points in concentration space that are reachable from a given initial condition, which is an affine subspace $$u_0 + S$$

When we say that a network is multistationary, or has the capacity for multiple steady states, what that means is that there is *some* stoichiometric compatibility class for which this is true. This means that there may be only one or zero steady states in some other stoichiometric compatibility class. The SCC can be inferred from the initial concentrations of the species in a reaction network via the network's conservation laws.

The summary function for the package is `networksummary`.

```@docs
networksummary
```

## Uniqueness of Steady States

The functionality that decides uniqueness of steady states is in the function `hasuniqueequilibria(rn)`. The uniqueness of steady states checks proceeds as follows. The function terminates every time a network is conclusively determined to have unique steady states or multi-stationarity.

 1. **Deficiency one and deficiency zero theorems**: a reaction network that satisfies the conditions of these two theorems will be guaranteed to have a unique steady state for every stoichiometric compatibility class.[^1]
 2. **Deficiency one and higher deficiency algorithm**: Both of these algorithms construct inequality systems that, if satsifiable, would indicate that a reaction network has the ability to admit multiple steady states in a given stoichiometric compatibility class.[^2]
 3. **Concordance**: This is a graph property that guarantees the uniqueness of steady states, since it guarantees that the species-formation rate function is injective. As a result, there is only one solution to $\frac{dc}{dt} = 0$.[^2]

Future desiderata for this check include methods from **monotone systems theory**,**degree theory**, and **embedded multi-stationary networks** - see the [Roadmap] for the package.

Other functionality related to the existence and uniqueness of steady states:

```@docs
# CatalystNetworkAnalysis.haspositivesteadystates
# CatalystNetworkAnalysis.hasperiodicsolutions
CatalystNetworkAnalysis.hasuniquesteadystates
# Catalyst.satisfiesdeficiencyzero
# Catalyst.satisfiesdeficiencyone
# Catalyst.iscomplexbalanced
# Catalyst.isdetailedbalanced
CatalystNetworkAnalysis.mixedvolume
isconcordant
deficiencyonealgorithm
```

Additionally, one might check the functions `iscomplexbalanced()` and `isdetailedbalanced()` implemented in the main Catalyst package, which work for kinetic systems (reaction networks with an assignment of parameter values). Both of these conditions are sufficient to guarantee that the reaction network will have unique steady states, for a given choice of parameters, and that these steady states will be *asymptotically stable* within each stoichiometric compatibility class.[^1]

## Stability

Bistability is a very important property of reaction networks, since switchlike behavior in biological networks is realized through bistability. The network analysis package currently does not have any tooling for stability analysis, but please see [this section](https://docs.sciml.ai/Catalyst/stable/steady_state_functionality/steady_state_stability_computation/) of the main Catalyst documentation for some tips.

## Persistence

A chemical reaction network is persistent if each species does not "die out" - that is, if the initial condition of the network has positive concentration in a given species, that species' concentration will never go to zero.

It is not currently possible to conclusively determine persistence in every case. The approach used by this package checks the necessary conditions and sufficient conditions outlined in this paper.[^4] Many persistent networks will be detectably persistent via this algorithm, and many non-persistent networks will be detectably non-persistent, but the function will return inconclusive in many cases.

```@docs
ispersistent
minimalsiphons
```

The Petri net algorithm relies on detecting siphons in the network, which follows a SAT-solver approach[^5] [^6].

### Concentration Robustness

A species in a reaction network has concentration robustness if its concentration in any steady state of the reaction network is the same. Species with this property will have the same concentration regardless of initial conditions or perturbations to the system.

 1. **Deficiency results**: these are very limited to deficiency one networks, in which robust species can be detected from the complex composition matrix. Additionally, there is a necessary deficiency-related condition for networks to have ACR, which can be used to detect networks without ACR[^3].
 2. **Algebraic results**: this class of results relies on results from real algebraic geometry. A network is concentration-robust if the locus of steady states lies entirely within a hyperplane where some species concentration is constant. Often, this boils down to finding basis elements of the form $(x_i - \alpha)$ for the steady state ideal[^8].

Like persistence, absolute concentration robustness currently cannot be detected conclusively in every case.

```@docs
isconcentrationrobust
```

Misc references[^9] [^10] [^11] [^12] [^13] [^14] [^15] [^16] [^17] [^18].

## References

[^1]: Feinberg, M. Foundations of Chemical Reaction Network Theory; Applied Mathematical Sciences; Springer International Publishing: Cham, 2019; Vol. 202. https://doi.org/10.1007/978-3-030-03858-8.
[^2]: Ji, H. Uniqueness of Equilibria for Complex Chemical Reaction Networks. Ph.D. Thesis, The Ohio State University, Columbus, OH, 2011.
[^3]: Eloundou-Mbebi, J. M. O.; Küken, A.; Omranian, N.; Kleessen, S.; Neigenfind, J.; Basler, G.; Nikoloski, Z. A Network Property Necessary for Concentration Robustness. Nat Commun 2016, 7 (1), 13255. https://doi.org/10.1038/ncomms13255.
[^4]: Angeli, D.; De Leenheer, P.; Sontag, E. D. A Petri Net Approach to the Study of Persistence in Chemical Reaction Networks. Math Biosci 2007, 210 (2), 598–618. https://doi.org/10.1016/j.mbs.2007.07.003.
[^5]: Shiu, A.; Sturmfels, B. Siphons in Chemical Reaction Networks. Bull. Math. Biol. 2010, 72 (6), 1448–1463. https://doi.org/10.1007/s11538-010-9502-y.
[^6]: Nabli, F.; Martinez, T.; Fages, F.; Soliman, S. On Enumerating Minimal Siphons in Petri Nets Using CLP and SAT Solvers: Theoretical and Practical Complexity. Constraints 2016, 21 (2), 251–276. https://doi.org/10.1007/s10601-015-9190-1.
[^8]: Puente, L. D. G.; Gross, E.; Harrington, H. A.; Johnston, M.; Meshkat, N.; Millán, M. P.; Shiu, A. Absolute Concentration Robustness: Algebra and Geometry. arXiv December 29, 2023. http://arxiv.org/abs/2401.00078 (accessed 2024-09-27).
[^9]: Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).
[^10]: Craciun, G.; Jin, J.; Yu, P. Y. An Efficient Characterization of Complex-Balanced, Detailed-Balanced, and Weakly Reversible Systems. SIAM J. Appl. Math. 2020, 80 (1), 183–205. https://doi.org/10.1137/19M1244494.
[^11]: Loman, T. E.; Ma, Y.; Ilin, V.; Gowda, S.; Korsbo, N.; Yewale, N.; Rackauckas, C.; Isaacson, S. A. Catalyst: Fast and Flexible Modeling of Reaction Networks.
[^12]: Angeli, D.; Ferrell, J. E.; Sontag, E. D. Detection of Multistability, Bifurcations, and Hysteresis in a Large Class of Biological Positive-Feedback Systems. Proceedings of the National Academy of Sciences 2004, 101 (7), 1822–1827. https://doi.org/10.1073/pnas.0308265100.
[^13]: Müller, S.; Regensburger, G. Generalized Mass Action Systems: Complex Balancing Equilibria and Sign Vectors of the Stoichiometric and Kinetic-Order Subspaces. SIAM J. Appl. Math. 2012, 72 (6), 1926–1947. https://doi.org/10.1137/110847056.
[^14]: Conradi, C.; Feliu, E.; Mincheva, M.; Wiuf, C. Identifying Parameter Regions for Multistationarity. PLoS Comput Biol 2017, 13 (10), e1005751. https://doi.org/10.1371/journal.pcbi.1005751.
[^15]: Craciun, G.; Nazarov, F.; Pantea, C. Persistence and Permanence of Mass-Action and Power-Law Dynamical Systems. SIAM J. Appl. Math. 2013, 73 (1), 305–329. https://doi.org/10.1137/100812355.
[^16]: Knight, D.; Shinar, G.; Feinberg, M. Sharper Graph-Theoretical Conditions for the Stabilization of Complex Reaction Networks. Math Biosci 2015, 262, 10–27. https://doi.org/10.1016/j.mbs.2015.01.002.
[^17]: Müller, S.; Feliu, E.; Regensburger, G.; Conradi, C.; Shiu, A.; Dickenstein, A. Sign Conditions for Injectivity of Generalized Polynomial Maps with Applications to Chemical Reaction Networks and Real Algebraic Geometry. arXiv October 30, 2014. http://arxiv.org/abs/1311.5493 (accessed 2024-08-15).
[^18]: Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.
