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

#include "python/python_float.h"

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
template<class T> void export_array();

void export_integer();
void export_rational();
template<class R> void export_float();
template<class R> void export_interval();
template<class R> void export_differential();
template<class R> void export_derivative();

template<class R> void export_vector();
template<class R> void export_matrix();
template<class R> void export_tensor();
template<class R> void export_linear_program();
template<class R> void export_interval_vector();
template<class R> void export_interval_matrix();
template<class R> void export_interval_tensor();

void export_binary_tree();
void export_lattice_set();
void export_lattice_map();

template<class R> void export_function();

template<class R> void export_point();
template<class R> void export_interval_point();
template<class R> void export_point_list();
template<class R> void export_set();
template<class R> void export_rectangle();
template<class R> void export_parallelotope();
template<class R> void export_simplex();
template<class R> void export_zonotope();
template<class R> void export_interval_zonotope();
template<class R> void export_polytope();
template<class R> void export_polyhedron();
template<class R> void export_empty_set();
template<class R> void export_rectangular_set();
template<class R> void export_polyhedral_set();
template<class R> void export_list_set();
template<class R> void export_grid();
template<class R> void export_grid_set();
template<class R> void export_partition_tree_set();
template<class R> void export_hybrid_set();

template<class R> void export_constraint();

template<class R> void export_map();
template<class R> void export_affine_map();
template<class R> void export_affine_multimap();
template<class R> void export_polynomial_map();
template<class R> void export_vector_field();
template<class R> void export_affine_vector_field();
template<class R> void export_set_based_hybrid_automaton();
template<class R> void export_constraint_based_hybrid_automaton();

template<class R> void export_evolution_parameters();
template<class R> void export_solve();
template<class R> void export_apply();
template<class R> void export_integrate();
template<class R> void export_detector();
template<class R> void export_set_based_hybrid_evolver();
template<class R> void export_constraint_based_hybrid_evolver();

void export_postscript_output();
void export_text_output();
void export_tex_output();

template<class R> void export_henon_map();
template<class R> void export_duffing_equation();
template<class R> void export_van_der_pol_equation();
template<class R> void export_lorenz_system();

using namespace Ariadne;
using namespace Ariadne::Base;
using namespace Ariadne::Python;

using Numeric::Integer;
using Numeric::Rational;

BOOST_PYTHON_MODULE(ariadne)
{
  export_debug();
  export_logging();
  export_exceptions();

  export_tribool();
  export_array<bool>();
  export_array<index_type>();
  export_array<size_type>();
  export_array<Integer>();
  export_array<Rational>();
  export_array<Float>();

  export_integer();
  export_rational();
  export_float<Float>();
  export_interval<Float>();
  export_differential<Rational>();
  export_differential<Float>();
  export_derivative<Rational>();
  export_derivative<Float>();

  export_vector<Float>();
  export_vector<Rational>();
  export_matrix<Float>();
  export_matrix<Rational>();
  export_tensor<Float>();
  export_tensor<Rational>();
  export_linear_program<Rational>();
  export_interval_vector<Float>();
  export_interval_matrix<Float>();
  export_interval_tensor<Float>();

  export_binary_tree();
  export_lattice_set();
  export_lattice_map();

  export_function<Rational>();
  export_function<Float>();

  export_set<Float>();
  export_point<Float>();
  export_interval_point<Float>();
  export_point_list<Float>();
  export_rectangle<Float>();
  export_zonotope<Float>();
  export_interval_zonotope<Float>();
  export_parallelotope<Float>();
  export_simplex<Float>();
  export_polytope<Float>();
  export_polyhedron<Float>();
  export_empty_set<Float>();
  export_rectangular_set<Float>();
  export_polyhedral_set<Float>();
  export_list_set<Float>();
  export_grid<Float>();
  export_grid_set<Float>();
  export_partition_tree_set<Float>();
  export_hybrid_set<Float>();

  export_constraint<Float>();

  export_map<Float>();
  export_affine_map<Float>();
  export_affine_multimap<Float>();
  export_polynomial_map<Float>();
  export_vector_field<Float>();
  export_affine_vector_field<Float>();
  export_set_based_hybrid_automaton<Float>();
  export_constraint_based_hybrid_automaton<Float>();

  export_evolution_parameters<Float>();
  export_solve<Float>();
  export_apply<Float>();
  export_integrate<Float>();
  export_detector<Float>();
  export_set_based_hybrid_evolver<Float>();
  export_constraint_based_hybrid_evolver<Float>();

  export_postscript_output();
  export_tex_output();
  export_text_output();

  export_henon_map<Float>();
  export_duffing_equation<Float>();
  export_van_der_pol_equation<Float>();
  export_lorenz_system<Float>();

}
