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

#include "real_typedef.h"

namespace Ariadne { namespace Numeric {
class Float64;
class MPFloat;
class Rational;
}}
using namespace Ariadne::Numeric;

#include "declarations.h"

#include <boost/python.hpp>

void export_tribool();
template<class T> void export_array();

void export_numeric();
template<class R> void export_function();
template<class R> void export_interval();

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

template<class R> void export_point();
template<class R> void export_interval_point();
template<class R> void export_point_list();
template<class R> void export_set();
template<class R> void export_rectangle();
template<class R> void export_parallelotope();
template<class R> void export_simplex();
template<class R> void export_zonotope();
template<class R> void export_polytope();
template<class R> void export_polyhedron();
template<class R> void export_polyhedral_set();
template<class R> void export_list_set();
template<class R> void export_grid_set();
template<class R> void export_partition_tree_set();
template<class R> void export_hybrid_set();

template<class R> void export_map();
template<class R> void export_affine_map();
template<class R> void export_affine_multimap();
template<class R> void export_polynomial_map();
template<class R> void export_vector_field();
template<class R> void export_affine_vector_field();
template<class R> void export_hybrid_automaton();

template<class R> void export_solve();
template<class R> void export_apply();
template<class R> void export_integrate();
template<class R> void export_hybrid_evolver();

void export_postscript_output();

template<class R> void export_henon_map();
template<class R> void export_duffing_equation();
template<class R> void export_van_der_pol_equation();
template<class R> void export_lorenz_system();

BOOST_PYTHON_MODULE(ariadne)
{
  export_tribool();
  export_array<bool>();
  export_array<Ariadne::index_type>();
  export_array<Ariadne::size_type>();
  export_array<Ariadne::Integer>();
  export_array<Ariadne::Rational>();
  export_array<Ariadne::Real>();

  export_numeric();
  export_function<MPFloat>();
  export_interval<Float64>();
  export_interval<MPFloat>();

  export_vector<Float64>();
  export_vector<MPFloat>();
  export_vector<Rational>();
  export_matrix<Float64>();
  export_matrix<MPFloat>();
  export_matrix<Rational>();
  export_tensor<Float64>();
  export_tensor<MPFloat>();
  export_tensor<Rational>();
  export_linear_program<Rational>();
  export_interval_vector<Float64>();
  export_interval_vector<MPFloat>();
  export_interval_matrix<Float64>();
  export_interval_matrix<MPFloat>();
  export_interval_tensor<Float64>();
  export_interval_tensor<MPFloat>();

  export_binary_tree();
  export_lattice_set();
  export_lattice_map();

  export_set<Ariadne::Real>();
  export_point<Ariadne::Real>();
  export_interval_point<Ariadne::Real>();
  export_point_list<Ariadne::Real>();
  export_rectangle<Ariadne::Real>();
  export_zonotope<Ariadne::Real>();
  export_parallelotope<Ariadne::Real>();
  export_simplex<Ariadne::Real>();
  export_polytope<Ariadne::Real>();
  export_polyhedron<Ariadne::Real>();
  export_polyhedral_set<Ariadne::Real>();
  export_list_set<Ariadne::Real>();
  export_grid_set<Ariadne::Real>();
  export_partition_tree_set<Ariadne::Real>();
  export_hybrid_set<Ariadne::Real>();

  export_map<Ariadne::Real>();
  export_affine_map<Ariadne::Real>();
  export_affine_multimap<Ariadne::Real>();
  export_polynomial_map<Ariadne::Real>();
  export_vector_field<Ariadne::Real>();
  export_affine_vector_field<Ariadne::Real>();
  export_hybrid_automaton<Ariadne::Real>();

  export_solve<Ariadne::Real>();
  export_apply<Ariadne::Real>();
  export_integrate<Ariadne::Real>();
  export_hybrid_evolver<Ariadne::Real>();

  export_postscript_output();

  export_henon_map<Ariadne::Real>();
  export_duffing_equation<Ariadne::Real>();
  export_van_der_pol_equation<Ariadne::Real>();
  export_lorenz_system<Ariadne::Real>();

}
