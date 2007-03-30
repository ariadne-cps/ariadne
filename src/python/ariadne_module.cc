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

#include "declarations.h"

#include <boost/python.hpp>

void export_debug();

void export_tribool();
template<class T> void export_array();

void export_integer();
void export_rational();
template<class R> void export_float();
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
template<class R> void export_interval_zonotope();
template<class R> void export_polytope();
template<class R> void export_polyhedron();
template<class R> void export_polyhedral_set();
template<class R> void export_list_set();
template<class R> void export_grid();
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
  export_debug();
  export_tribool();
  export_array<bool>();
  export_array<Ariadne::index_type>();
  export_array<Ariadne::size_type>();
  export_array<Ariadne::Integer>();
  export_array<Ariadne::Rational>();
  export_array<Ariadne::Float>();

  export_integer();
  export_rational();
  export_float<Ariadne::Float>();
  export_interval<Ariadne::Float>();

  export_vector<Ariadne::Float>();
  export_vector<Ariadne::Rational>();
  export_matrix<Ariadne::Float>();
  export_matrix<Ariadne::Rational>();
  export_tensor<Ariadne::Float>();
  export_tensor<Ariadne::Rational>();
  export_linear_program<Ariadne::Rational>();
  export_interval_vector<Ariadne::Float>();
  export_interval_matrix<Ariadne::Float>();
  export_interval_tensor<Ariadne::Float>();

  export_binary_tree();
  export_lattice_set();
  export_lattice_map();

  export_set<Ariadne::Float>();
  export_point<Ariadne::Float>();
  export_interval_point<Ariadne::Float>();
  export_point_list<Ariadne::Float>();
  export_rectangle<Ariadne::Float>();
  export_zonotope<Ariadne::Float>();
  export_interval_zonotope<Ariadne::Float>();
  export_parallelotope<Ariadne::Float>();
  export_simplex<Ariadne::Float>();
  export_polytope<Ariadne::Float>();
  export_polyhedron<Ariadne::Float>();
  export_polyhedral_set<Ariadne::Float>();
  export_list_set<Ariadne::Float>();
  export_grid<Ariadne::Float>();
  export_grid_set<Ariadne::Float>();
  export_partition_tree_set<Ariadne::Float>();
  export_hybrid_set<Ariadne::Float>();

  export_map<Ariadne::Float>();
  export_affine_map<Ariadne::Float>();
  export_affine_multimap<Ariadne::Float>();
  export_polynomial_map<Ariadne::Float>();
  export_vector_field<Ariadne::Float>();
  export_affine_vector_field<Ariadne::Float>();
  export_hybrid_automaton<Ariadne::Float>();

  export_solve<Ariadne::Float>();
  export_apply<Ariadne::Float>();
  export_integrate<Ariadne::Float>();
  export_hybrid_evolver<Ariadne::Float>();

  export_postscript_output();

  export_henon_map<Ariadne::Float>();
  export_duffing_equation<Ariadne::Float>();
  export_van_der_pol_equation<Ariadne::Float>();
  export_lorenz_system<Ariadne::Float>();

}
