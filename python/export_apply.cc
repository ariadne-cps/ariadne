/***************************************************************************
 *            python/export_apply.cc
 *
 *  6 February 2006
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

#include "apply.h"
#include "henon_map.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

#include "real_typedef.h"

using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;

inline Parallelopiped<Real> apply_to_parallelopiped(const Map<Real>& f, const Parallelopiped<Real>& p) {
  return apply(f,p); 
}

inline ListSet<Real,Parallelopiped> apply_to_parallelopiped_list_set(const Map<Real>& f, const ListSet<Real,Parallelopiped>& pls) {
  return apply(f,pls); 
}

inline GridMaskSet<Real> 
chainreach_of_rectangle_list_set(const Map<Real>& f, 
                                 const ListSet<Real,Rectangle>& rls,
                                 const FiniteGrid<Real>& g,
                                 const Rectangle<Real>& bb)
{
  return chainreach(f,rls,g,bb); 
}

inline ListSet<Real,Parallelopiped> apply_henon_map_to_parallelopiped_list_set(const HenonMap<Real>& f, const ListSet<Real,Parallelopiped>& pls) {
  std::cerr << "apply_henon_map_to_parallelopiped_list_set" << std::endl << std::flush;
  return apply(f,pls); 
}


void export_apply() {
  def("apply", apply_to_parallelopiped, "apply the image of a polynomial to a set" );
  def("apply", apply_to_parallelopiped_list_set, "apply the image of a polynomial to a set" );
  def("chainreach", &chainreach_of_rectangle_list_set, "chain reach of a set" );
//  def("apply", apply_henon_map_to_parallelopiped_list_set, "apply the image of a polynomial to a set" );
}
