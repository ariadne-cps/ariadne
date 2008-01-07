/***************************************************************************
 *            python/export_approximator.cc
 *
 *  Copyright  2007  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is diself_ns::stributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Library General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "python/float.h"

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "system/map.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/standard_approximator.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS>
class ApproximatorWrapper
  : public ApproximatorInterface<BS>,
    public wrapper< ApproximatorInterface<BS> >
{
  typedef typename BS::real_type R;
  typedef Interval<R> I;
 public:
  ApproximatorWrapper<BS>* clone() const { return this->get_override("clone")(); }
  BS basic_set(const Box<R>&) const { return this->get_override("basic_set")(); }
  R radius(const BS&) const { return this->get_override("radius")(); }
  Box<R> bounding_box(const BS&) const { return this->get_override("bounding_box")(); }
  GridCellListSet<R> outer_approximation(const BS&, const Grid<R>&) const { return this->get_override("outer_approximation")(); }
};


template<class R>
void export_approximator() 
{
  typedef Zonotope<R> ZBS;

  class_< ApproximatorWrapper<ZBS>, boost::noncopyable >("ZonotopeApproximatorInterface",init<>());

  class_< StandardApproximator<ZBS>, bases< ApproximatorInterface<ZBS> > >
    approximator_class("StandardApproximator",init<>());

}

template void export_approximator<FloatPy>();
