# Testing of the linear programming utilities in CatalystNetworkAnalysis

import CatalystNetworkAnalysis as C
# Test that matrices with positive element in their image space are identified.
let
    a = rand(1:10, 10, 10)
    @test has_positive_solution(a)
    
    b = -rand(1:10, 10, 10)
    @test has_positive_solution(b)

    a = [1 -1; -1 1]
    @test !has_positive_solution(b)

    b = [-1 1 1 1;
         1 -1 1 1;
         1 1 -1 1
         1 1 1 -1]
    @test !has_positive_solution(b)
    @test has_positive_solution(b, nonneg=true)
end

# Testing that extreme rays are properly identified.
let
end

# Testing elementary flux modes are found.
let
end
