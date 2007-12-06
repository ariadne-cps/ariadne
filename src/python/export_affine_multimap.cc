/***************************************************************************
 *            python/export_affine_map_with_set.cc
 *
 *  31 May 2006
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/simplex.h"
#include "geometry/zonotope.h"
#include "geometry/polytope.h"
#include "system/affine_multimap.h"
#include "system/affine_map.h"
#include "python/utilities.h"

using namespace Ariadne;
using namespace Ariadne::LinearAlgebra;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class R>
void export_affine_multimap() 
{
  typedef AffineMultiMap<R,Geometry::Rectangle> RAffineMultiMapRectangle;
  typedef Matrix<R> RMatrix;
  typedef Geometry::Rectangle<R> RRectangle;
  
  class_< AffineMultiMap<R,Rectangle> >("AffineMultiMapRectangle",init< Matrix<R>, Rectangle<R> >())
    .def(init< AffineMultiMap<R,Rectangle> >())
    .def("argument_dimension",&AffineMultiMap<R,Rectangle>::argument_dimension)
    .def("result_dimension",&AffineMultiMap<R,Rectangle>::result_dimension)
//    .def("__call__",(Rectangle<R>(AffineMultiMap<R,Rectangle>::*)(const Point<R>&)const)(AffineMultiMap<R,Rectangle>::operator()))
//    .def("__call__",(Rectangle<R>(AffineMultiMap<R,Rectangle>::*)(const Rectangle<R>&)const)(AffineMultiMap<R,Rectangle>::operator()))
  ;

  class_< AffineMultiMap<R,Zonotope> >("AffineMultiMapZonotope",init< Matrix<R>, Zonotope<R> >())
    .def(init< AffineMultiMap<R,Zonotope> >())
    .def("argument_dimension",&AffineMultiMap<R,Zonotope>::argument_dimension)
    .def("result_dimension",&AffineMultiMap<R,Zonotope>::result_dimension)
    //.def("__call__",(Zonotope<R>(AffineMultiMap<R,Zonotope>::*)(const Point<R>&)const)(AffineMultiMap<R,Zonotope>::operator()))
    //.def("__call__",(Zonotope<R>(AffineMultiMap<R,Zonotope>::*)(const Zonotope<R>&)const)(AffineMultiMap<R,Zonotope>::operator()))
  ;

}

template void export_affine_multimap<FloatPy>();
