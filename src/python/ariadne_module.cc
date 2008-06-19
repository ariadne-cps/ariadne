/***************************************************************************
 *            python/ariadne_module.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "base/types.h"
#include "base/declarations.h"
#include "numeric/declarations.h"
#include "linear_algebra/declarations.h"
#include "combinatoric/declarations.h"
#include "geometry/declarations.h"
#include "system/declarations.h"
#include "evaluation/declarations.h"

#include <boost/python.hpp>

void export_debug();
void export_logging();
void export_exceptions();

void export_tribool();
void export_arrays();

void export_integer();
void export_rational();
template<class R> void export_float();
template<class R> void export_interval();

template<class R> void export_vector();
template<class R> void export_covector();
template<class R> void export_matrix();
template<class R> void export_linear_program();
template<class R> void export_interval_vector();
template<class R> void export_interval_covector();
template<class R> void export_interval_matrix();

void export_multi_index();
void export_binary_tree();
void export_lattice_set();
void export_lattice_map();
void export_subdivision_set();

template<class R> void export_function_interface();
template<class R> void export_python_function();
template<class R> void export_affine_function();
template<class R> void export_affine_variable();
template<class R> void export_affine_derivative();
template<class R> void export_taylor_series();
template<class R> void export_taylor_variable();
template<class R> void export_taylor_derivative();
template<class R> void export_flow_model();
template<class R> void export_taylor_model();
template<class R> void export_polynomial_variable();

template<class R> void export_point();
template<class R> void export_interval_point();
template<class R> void export_point_list();
template<class R> void export_segment();
template<class R> void export_interpolated_curve();
template<class R> void export_box();
template<class R> void export_box_list_set();
template<class R> void export_set();
template<class R> void export_interval_set();
template<class R> void export_rectangle();
template<class R> void export_simplex();
template<class R> void export_zonotope();
template<class R> void export_polytope();
template<class R> void export_polyhedron();
template<class R> void export_empty_set();
template<class R> void export_list_set();
template<class R> void export_grid();
template<class R> void export_grid_set();
template<class R> void export_partition_tree_set();
template<class R> void export_constraint_set();
template<class R> void export_rectangular_set();
template<class R> void export_polyhedral_set();
template<class R> void export_hybrid_set();

template<class R> void export_constraint();

template<class R> void export_orbit();

template<class R> void export_approximator();
template<class R> void export_subdivider();
template<class R> void export_reducer();

template<class R> void export_map();
template<class R> void export_vector_field();
template<class R> void export_hybrid_automaton();
template<class R> void export_constraint_based_hybrid_automaton();

template<class R> void export_evolution_parameters();
template<class R> void export_solver();
template<class R> void export_applicator();
template<class R> void export_integrator();
template<class R> void export_detector();
template<class R> void export_discretiser();
template<class R> void export_evolver();
template<class R> void export_reachability_analyser();
template<class R> void export_model_checker();

void export_text_output();
void export_latex_output();
void export_postscript_output();

template<class R> void export_henon_map();
template<class R> void export_duffing_equation();
template<class R> void export_van_der_pol_equation();
template<class R> void export_lorenz_system();

using namespace Ariadne;
using namespace Ariadne::Python;

BOOST_PYTHON_MODULE(ariadne)
{
  export_debug();
  export_logging();
  export_exceptions();

  export_tribool();
  export_arrays();

  export_integer();
  export_rational();
  export_float<FloatPy>();
  export_interval<FloatPy>();

  export_vector<FloatPy>();
  export_vector<Rational>();
  export_covector<FloatPy>();
  export_covector<Rational>();
  export_matrix<FloatPy>();
  export_matrix<Rational>();
  export_interval_vector<FloatPy>();
  export_interval_covector<FloatPy>();
  export_interval_matrix<FloatPy>();

  export_linear_program<Rational>();

  export_multi_index();
  export_binary_tree();
  export_lattice_set();
  export_lattice_map();
  export_subdivision_set();

  export_function_interface<Rational>();
  export_function_interface<FloatPy>();
  export_python_function<Rational>();
  export_python_function<FloatPy>();
  export_affine_function<Rational>();
  export_affine_function<FloatPy>();
  //export_affine_variable<Rational>();
  //export_affine_variable<FloatPy>();
  //export_affine_derivative<Rational>();
  //export_affine_derivative<FloatPy>();
  export_taylor_series<FloatPy>();
  export_taylor_variable<FloatPy>();
  export_taylor_derivative<FloatPy>();
  export_taylor_model<FloatPy>();
  //export_polynomial_variable<FloatPy>();
  export_flow_model<FloatPy>();
  //export_taylor_series<Rational>();
  //export_taylor_variable<Rational>();
  //export_taylor_derivative<Rational>();

  export_point<Rational>();
  export_segment<Rational>();
  export_interpolated_curve<Rational>();
  export_polytope<Rational>();
  export_polyhedron<Rational>();

  export_set<FloatPy>();
  export_point<FloatPy>();
  export_interval_point<FloatPy>();
  export_point_list<FloatPy>();
  export_segment<FloatPy>();
  export_interpolated_curve<FloatPy>();
  export_interval_set<FloatPy>();
  export_box<FloatPy>();
  export_box_list_set<FloatPy>();
  export_rectangle<FloatPy>();
  export_zonotope<FloatPy>();
  export_simplex<FloatPy>();
  export_polytope<FloatPy>();
  export_polyhedron<FloatPy>();
  export_empty_set<FloatPy>();
  export_constraint_set<FloatPy>();
  export_rectangular_set<FloatPy>();
  export_polyhedral_set<FloatPy>();
  export_list_set<FloatPy>();
  export_grid<FloatPy>();
  export_grid_set<FloatPy>();
  export_partition_tree_set<FloatPy>();
  export_hybrid_set<FloatPy>();

  export_constraint<FloatPy>();

  export_orbit<FloatPy>();

  export_approximator<FloatPy>();
  export_subdivider<FloatPy>();
  export_reducer<FloatPy>();

  export_map<FloatPy>();
  export_vector_field<FloatPy>();
  export_hybrid_automaton<FloatPy>();
  //export_constraint_based_hybrid_automaton<FloatPy>();

  export_evolution_parameters<FloatPy>();
  export_solver<FloatPy>();
  export_applicator<FloatPy>();
  export_integrator<FloatPy>();
  export_evolver<FloatPy>();
  export_reachability_analyser<FloatPy>();
  export_model_checker<FloatPy>();

  export_text_output();
  export_latex_output();
  export_postscript_output();

  export_henon_map<FloatPy>();
  export_duffing_equation<FloatPy>();
  export_van_der_pol_equation<FloatPy>();
  export_lorenz_system<FloatPy>();

}
