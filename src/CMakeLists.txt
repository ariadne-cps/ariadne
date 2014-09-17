message(${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
# Non-portable fix to avoid user header directories being searched for by system headers
add_definitions(-iquote${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(Boost 1.36.0)

add_library(ariadne
    operators.cc
    integer.cc
    rational.cc
    decimal.cc
    float.cc
    float-approximate.cc
    float-validated.cc
    interval.cc
    real.cc
    vector.cc
    matrix.cc
    differential.cc
    affine.cc
    polynomial.cc
    expression.cc
    function.cc
    algebra.cc
    affine_model.cc
    taylor_model.cc
    taylor_function.cc
    paver.cc
    solver.cc
    integrator.cc
    constraint_solver.cc
    simplex_algorithm.cc
    linear_programming.cc
    nonlinear_programming.cc
    point.cc
    box.cc
    curve.cc
    orbit.cc
    function_set.cc
    grid_set.cc
    affine_set.cc
    enclosure.cc
    expression_set.cc
    map_evolver.cc
    paver.cc
    vector_field_evolver.cc
    logging.cc
    graphics.cc
    drawer.cc
    hybrid_set.cc
    hybrid_enclosure.cc
    hybrid_automaton-monolithic.cc
    hybrid_automaton-composite.cc
    hybrid_automaton-restrictive.cc
    hybrid_simulator.cc
    hybrid_evolver.cc
    hybrid_reachability_analyser.cc
    hybrid_graphics.cc
)

add_library(ariadne-extra
    runge_kutta_integrator.cc
    geometry.cc
    zonotope.cc
    polytope.cc
    polyhedron.cc
    textplot.cc
)

target_link_libraries(ariadne gmpxx gmp)
#FIXME: Make boost_system non-required
target_link_libraries(ariadne boost_system)
#FIXME: Make Cairo and Boost Serialization optional dependencies
target_link_libraries(ariadne cairo)
target_link_libraries (ariadne boost_serialization)
target_link_libraries(ariadne-extra ariadne)

