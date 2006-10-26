/***************************************************************************
 *            python/geometry_module.cc
 *
 *  21 October 2005
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

#include <boost/python.hpp>

template<class R> void export_point();
template<class R> void export_point_list();
template<class R> void export_rectangle();
template<class R> void export_parallelotope();
template<class R> void export_simplex();
template<class R> void export_zonotope();
template<class R> void export_polytope();
template<class R> void export_polyhedron();
template<class R> void export_list_set();
template<class R> void export_grid_set();
template<class R> void export_partition_tree_set();

BOOST_PYTHON_MODULE(geometry)
{
  export_point<Ariadne::Real>();
  export_point_list<Ariadne::Real>();
  export_rectangle<Ariadne::Real>();
  export_zonotope<Ariadne::Real>();
  export_parallelotope<Ariadne::Real>();
  export_simplex<Ariadne::Real>();
  export_polytope<Ariadne::Real>();
  export_polyhedron<Ariadne::Real>();
  export_list_set<Ariadne::Real>();
  export_grid_set<Ariadne::Real>();
  export_partition_tree_set<Ariadne::Real>();
}
