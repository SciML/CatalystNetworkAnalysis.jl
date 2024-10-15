# This file contains functionality for translated reaction networks, which extend the results
# of deficiency theory to a much wider range of networks. 

# (1) Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.
 
# Functionality for translating chemical reaction networks. 

"""
    translate(rn::ReactionNetwork)

    Given a reaction network, attempt to find a strongly-resolvable translation that is weakly-reversible and deficiency zero. Such translations are useful because they can be easily parameterized, and their multistability characteristics can be easily determined. See  
"""
function WRDZtranslation(rn::ReactionNetwork) 
    # Get the complexes and define generalized complexes
    kineticcomplexes = complexstoichmat(rn)
    translatedcomplexes = copy(complexstoichmat(rn))

    # Find elementary flux modes and reactions with shared source complexes
    EFMs = elementaryfluxmodes(rn)
    # Check that EFMs are unitary and cover R

    # Solve BLP
    model = Model(HiGHS.Optimizer);
    set_silent(model)
    set_optimizer_attribute(model, "mip_feasibility_tolerance", 1e-10)
    @variable(model, edge[1:r, 1:r], Bin)
    @objective(model, Min, sum(edge))
    # Add CS, EM, PS constraints
end


# (1) Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).

# Generating the parameterization from the above paper. 

# Kinetic deficiency zero reaction networks have positive parameterizations. 
function kineticdeficiency(rn::ReactionNetwork) 
end

"""
    Given a reaction system, compute the positive parameterization of the system. 
    1) For generalized mass action systems with kinetic deficiency zero. 
    2) Algorithm for positive kinetic deficiency systems. 
    3) Monomial parameterization for systems with toric steady states (Millan et al. 2012)
    4) Using Matroids (fall-back), checking subsets of d variables? 

    The output of this function is a symbolic function representing the parameterization. It takes some subset of the species and maps them into a steady state.  
"""
function positiveparameterization(rn::ReactionSystem; variable_search = false) 
    if kineticdeficiency(rn) == 0
    end

    # κ ∘ τ
    # Find a spanning forest
    # Find B such that im B = ker M^T
    # H is a generalized inverse of M
end

