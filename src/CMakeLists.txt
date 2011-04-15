message(${PROJECT_SOURCE_DIR})

include_directories(${PROJECT_SOURCE_DIR})
# Non-portable fix to avoid user header directories being searched for by system headers
add_definitions(-iquote${PROJECT_SOURCE_DIR}/include)
#include_directories(${PROJECT_SOURCE_DIR}/include)

find_package(Boost 1.36.0) 

add_library(ariadne 
    integer.cc 
	rational.cc 
	float.cc 
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
	taylor_model.cc  
	taylor_function.cc  
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
	taylor_set.cc 
	grid_set.cc 
	affine_set.cc 
	map_evolver.cc 
	vector_field_evolver.cc 
	discretiser.cc 
	logging.cc 
	graphics.cc 
	hybrid_set.cc 
	hybrid_enclosure.cc 
	hybrid_automaton-monolithic.cc 
	hybrid_automaton-composite.cc 
	hybrid_evolver.cc 
	hybrid_discretiser.cc 
	hybrid_reachability_analyser.cc
)

add_library(ariadne-extra
	runge_kutta_integrator.cc
	geometry.cc
	zonotope.cc
	polytope.cc
	polyhedron.cc
	textplot.cc
)

target_link_libraries(ariadne gmpxx)
#FIXME: Make Cairo and Boost Serialization optional dependencies
target_link_libraries(ariadne cairo)
target_link_libraries (ariadne boost_serialization)
target_link_libraries(ariadne-extra ariadne)

