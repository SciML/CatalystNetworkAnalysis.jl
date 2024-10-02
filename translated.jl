### This file contains functionality for translated reaction networks, which extend the results
# of deficiency theory to a much wider range of networks. 

# (1) Johnston, M. D. Translated Chemical Reaction Networks. Bull Math Biol 2014, 76 (5), 1081–1116. https://doi.org/10.1007/s11538-014-9947-5.
 
# Functionality for translating chemical reaction networks. 
function translate(rn::ReactionNetwork) 
    
end


# (1) Johnston, M. D.; Müller, S.; Pantea, C. A Deficiency-Based Approach to Parametrizing Positive Equilibria of Biochemical Reaction Systems. arXiv May 23, 2018. http://arxiv.org/abs/1805.09295 (accessed 2024-09-09).

# Generating the parameterization from the above paper. 

# Kinetic deficiency zero reaction networks have positive parameterizations. 
function kineticdeficiency(rn::ReactionNetwork) 
    
end

function positiveparameterization(rn::ReactionSystem) 
    
end
