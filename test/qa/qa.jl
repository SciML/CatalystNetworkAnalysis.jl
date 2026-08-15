using SciMLTesting, CatalystNetworkAnalysis, JET
using Test

# aqua_broken sub-checks (tracked in issue #70):
#   stale_deps: ReactionNetworkImporters, PolynomialRoots, ModelingToolkit, SBMLToolkit
#
# jet_broken: ~14 possible errors (report mode + @test_broken) — issue #70
#
# ExplicitImports:
#   no_implicit_imports is ei_broken — the package `using`s several large
#   dependencies (Catalyst, Oscar, JuMP, Graphs, DynamicPolynomials,
#   Satisfiability, ...) and relies on ~75 of their exported names implicitly;
#   making every one explicit is a large refactor tracked in issue #70.
#   The qualified-access ignore-lists below are names that are not public in
#   their *source* packages (other people's packages) but are accessed
#   intentionally via Mod.name; ignored per ExplicitImports' public-API checks.
run_qa(
    CatalystNetworkAnalysis;
    aqua_broken = (:stale_deps,),
    jet_broken = true,
    ei_kwargs = (;
        all_qualified_accesses_are_public = (;
            ignore = (
                :Library,              # CDDLib
                :Optimizer,            # HiGHS
                :assemble_oderhs,      # Catalyst
                :cycles,               # Catalyst
                :get_networkproperties, # Catalyst
                :n_positive_roots,     # Hecke
                :symmap_to_varmap,     # Catalyst
            ),
        ),
    ),
    ei_broken = (:no_implicit_imports,),
)
