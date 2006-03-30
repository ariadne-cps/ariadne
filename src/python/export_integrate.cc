/***************************************************************************
 *            python/export_integrate.cc
 *
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "evaluation/integrate.h"

#include <boost/python.hpp>

using boost::python::class_;
using boost::python::init;
using boost::python::self;
using boost::python::def;
using boost::python::return_value_policy;
using boost::python::copy_const_reference;

#include "python/real_typedef.h"

using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;

typedef Ariadne::Interval<Real> RInterval;
typedef Ariadne::Geometry::Rectangle<Real> RRectangle;
typedef Ariadne::Geometry::Parallelotope<Real> RParallelotope;
typedef Ariadne::Geometry::ListSet<Real,Parallelotope> RParallelotopeListSet;
typedef Ariadne::Evaluation::VectorField<Real> RVectorField;


/*
inline GridMaskSet<Real> 
chainreach_of_rectangle_list_set(const Map<Real>& f, 
                                 const ListSet<Real,Rectangle>& rls,
                                 const FiniteGrid<Real>& g,
                                 const Rectangle<Real>& bb)
{
  return chainreach(f,rls,g,bb); 
}
*/


void export_integrate() {
  typedef RRectangle (*IntRectFunc) (const RVectorField&, const RRectangle&, const Real&);
  typedef RParallelotope (*IntPltpFunc) (const RVectorField&, const RParallelotope&, const Real&);
  typedef RParallelotope (*IntPltpIntvlFunc) (const RVectorField&, const RParallelotope&, const RInterval&);
  typedef RParallelotopeListSet (*IntLSPltpFunc) (const RVectorField&, const RParallelotopeListSet&, const Real&);
  
  def("integrate", IntRectFunc(&integrate), "integrate a vector field over a set");
  def("integrate", IntPltpFunc(&integrate));
  def("integrate", IntPltpIntvlFunc(&integrate));
  def("integrate_to", IntPltpFunc(&integrate_to));
  def("integrate", IntLSPltpFunc(&integrate));
//  def("chainreach", &chainreach_of_rectangle_list_set, "chain reach of a set" );
}
