### Network Analysis Algorithms Overview
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
Functions that require one to supply a set of parameters will take `p`, a dictionary, vector, or
tuple of symbol-value pairs as an argument. 

Another key notion is that of the **stoichiometric compatibility class**, which is the set of points 
in concentration space that are reachable from a given initial condition. When we say that a network
is multistationary, or has the capacity for multiple steady states, what that means is that there is
*some* stoichiometric compatibility class for which this is true. This means that there may be only one
or zero steady states in some other stoichiometric compatibility class. The SCC can be inferred from 
the initial concentrations of the species in a reaction network. Functions that optionally allow one
to specify an initial set of concentrations will take `u0`, a dictionary, vector, or tuple of symbol-value
pairs as an argument. 

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

**In progress**: 
1. **Monotone systems theory**: 
2. **Degree theory**: 
3. **Embedded multi-stationary networks**: 

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
A species in a reaction network has concentration robustness if its concentration in any 
steady state of the reaction network is the same. Species with this property will have
the same concentration regardless of initial conditions or perturbations to the system. 

1. Deficiency results: these are very limited to deficiency one networks, in which robust species can be detected from the complex composition matrix.
2. Algebraic results: this class of results relies on results from real algebraic geometry, in which concentration robustness is equivalent to finding 
basis elements of the form $(x_i - \alpha)$ for the steady state ideal.

Like persistence, absolute concentration robustness currently cannot be detected conclusively in every case. 

### References
(1) Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).
(2) Hernandez, B. S.; Lubenia, P. V. N.; Johnston, M. D.; Kim, J. K. A Framework for Deriving Analytic Steady States of Biochemical Reaction Networks. PLOS Computational Biology 2023, 19 (4), e1011039. https://doi.org/10.1371/journal.pcbi.1011039.
(3) Gillespie, D. T. A Rigorous Derivation of the Chemical Master Equation. Physica A: Statistical Mechanics and its Applications 1992, 188 (1–3), 404–425. https://doi.org/10.1016/0378-4371(92)90283-V.
(4) Puente, L. D. G.; Gross, E.; Harrington, H. A.; Johnston, M.; Meshkat, N.; Millán, M. P.; Shiu, A. Absolute Concentration Robustness: Algebra and Geometry. arXiv December 29, 2023. http://arxiv.org/abs/2401.00078 (accessed 2024-09-27).
(5) Gross, E.; Harrington, H. A.; Rosen, Z.; Sturmfels, B. Algebraic Systems Biology: A Case Study for the Wnt Pathway. arXiv February 11, 2015. http://arxiv.org/abs/1502.03188 (accessed 2024-10-11).
(6) Craciun, G.; Jin, J.; Yu, P. Y. An Efficient Characterization of Complex-Balanced, Detailed-Balanced, and Weakly Reversible Systems. SIAM J. Appl. Math. 2020, 80 (1), 183–205. https://doi.org/10.1137/19M1244494.
(7) Schnoerr, D.; Sanguinetti, G.; Grima, R. Approximation and Inference Methods for Stochastic Biochemical Kinetics—a Tutorial Review. J. Phys. A: Math. Theor. 2017, 50 (9), 093001. https://doi.org/10.1088/1751-8121/aa54d9.
(8) Ji, W.; Deng, S. Autonomous Discovery of Unknown Reaction Pathways from Data by Chemical Reaction Neural Network.
(9) Loman, T. E.; Ma, Y.; Ilin, V.; Gowda, S.; Korsbo, N.; Yewale, N.; Rackauckas, C.; Isaacson, S. A. Catalyst: Fast and Flexible Modeling of Reaction Networks.
(10) Rosen, Z. Computing Algebraic Matroids. arXiv April 8, 2014. http://arxiv.org/abs/1403.8148 (accessed 2024-10-11).
(11) Sanft, K. R.; Othmer, H. G. Constant-Complexity Stochastic Simulation Algorithm with Optimal Binning. The Journal of Chemical Physics 2015, 143 (7), 074108. https://doi.org/10.1063/1.4928635.
(12) Angeli, D.; Ferrell, J. E.; Sontag, E. D. Detection of Multistability, Bifurcations, and Hysteresis in a Large Class of Biological Positive-Feedback Systems. Proceedings of the National Academy of Sciences 2004, 101 (7), 1822–1827. https://doi.org/10.1073/pnas.0308265100.
(13) Zanghellini, J.; Ruckerbauer, D. E.; Hanscho, M.; Jungreuthmayer, C. Elementary Flux Modes in a Nutshell: Properties, Calculation and Applications. Biotechnology Journal 2013, 8 (9), 1009–1016. https://doi.org/10.1002/biot.201200269.
(14) Feinberg, M. Foundations of Chemical Reaction Network Theory; Applied Mathematical Sciences; Springer International Publishing: Cham, 2019; Vol. 202. https://doi.org/10.1007/978-3-030-03858-8.
(15) Carlson, R.; Srienc, F. Fundamental Escherichia Coli Biochemical Pathways for Biomass and Energy Production: Creation of Overall Flux States. Biotechnol Bioeng 2004, 86 (2), 149–162. https://doi.org/10.1002/bit.20044.
(16) Müller, S.; Regensburger, G. Generalized Mass Action Systems: Complex Balancing Equilibria and Sign Vectors of the Stoichiometric and Kinetic-Order Subspaces. SIAM J. Appl. Math. 2012, 72 (6), 1926–1947. https://doi.org/10.1137/110847056.
(17) Conradi, C.; Feliu, E.; Mincheva, M.; Wiuf, C. Identifying Parameter Regions for Multistationarity. PLoS Comput Biol 2017, 13 (10), e1005751. https://doi.org/10.1371/journal.pcbi.1005751.
(18) Bihan, F.; Dickenstein, A.; Giaroli, M. Lower Bounds for Positive Roots and Regions of Multistationarity in Chemical Reaction Networks. Journal of Algebra 2020, 542, 367–411. https://doi.org/10.1016/j.jalgebra.2019.10.002.
(19) Craciun, G.; Joshi, B.; Pantea, C.; Tan, I. Multistationarity in Cyclic Sequestration-Transmutation Networks. arXiv April 11, 2022. http://arxiv.org/abs/2110.13975 (accessed 2024-09-18).
(20) Huck, W. T. S. Oscillating Chemical Reaction Networks Stopped Cold. Nat Chem Eng 2024, 1–2. https://doi.org/10.1038/s44286-024-00092-8.
(21) Craciun, G.; Nazarov, F.; Pantea, C. Persistence and Permanence of Mass-Action and Power-Law Dynamical Systems. SIAM J. Appl. Math. 2013, 73 (1), 305–329. https://doi.org/10.1137/100812355.
(22) Kauffman, S.; Peterson, C.; Samuelsson, B.; Troein, C. Random Boolean Network Models and the Yeast Transcriptional Network. Proceedings of the National Academy of Sciences 2003, 100 (25), 14796–14799. https://doi.org/10.1073/pnas.2036429100.
(23) Ferrell, J. E. Self-Perpetuating States in Signal Transduction: Positive Feedback, Double-Negative Feedback and Bistability. Curr Opin Cell Biol 2002, 14 (2), 140–148. https://doi.org/10.1016/s0955-0674(02)00314-9.
(24) Knight, D.; Shinar, G.; Feinberg, M. Sharper Graph-Theoretical Conditions for the Stabilization of Complex Reaction Networks. Math Biosci 2015, 262, 10–27. https://doi.org/10.1016/j.mbs.2015.01.002.
(25) Müller, S.; Feliu, E.; Regensburger, G.; Conradi, C.; Shiu, A.; Dickenstein, A. Sign Conditions for Injectivity of Generalized Polynomial Maps with Applications to Chemical Reaction Networks and Real Algebraic Geometry. arXiv October 30, 2014. http://arxiv.org/abs/1311.5493 (accessed 2024-08-15).
(26) Banaji, M. Some Bounds on Positive Equilibria in Mass Action Networks. arXiv September 10, 2024. http://arxiv.org/abs/2409.06877 (accessed 2024-09-18).
(27) Anderson, D. F.; Cappelletti, D.; Kim, J.; Nguyen, T. D. Tier Structure of Strongly Endotactic Reaction Networks. Stochastic Processes and their Applications 2020, 130 (12), 7218–7259. https://doi.org/10.1016/j.spa.2020.07.012.
(28) Telek, M. L.; Feliu, E. Topological Descriptors of the Parameter Region of Multistationarity: Deciding upon Connectivity. PLOS Computational Biology 2023, 19 (3), e1010970. https://doi.org/10.1371/journal.pcbi.1010970.
(29) Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.




