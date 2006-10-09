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

#include <boost/python.hpp>
#include "numeric/float64.h"
#include "numeric/mpfloat.h"
#include "numeric/rational.h"


template<typename R> void export_point();
template<typename R> void export_point_list();
template<typename R> void export_rectangle();
void export_parallelotope();
void export_simplex();
void export_zonotope();
void export_polytope();
void export_polyhedron();
void export_list_set();
void export_grid_set();
void export_partition_tree_set();

BOOST_PYTHON_MODULE(geometry)
{
  export_point<Ariadne::Numeric::MPFloat>();
  export_point_list<Ariadne::Numeric::MPFloat>();
  export_rectangle<Ariadne::Numeric::MPFloat>();
  export_zonotope();
  export_parallelotope();
  export_simplex();
  export_polytope();
  export_polyhedron();
  export_list_set();
  export_grid_set();
  export_partition_tree_set();
}
