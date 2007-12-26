/***************************************************************************
 *            python/export_applicator.cc
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
#include "system/map_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/standard_applicator.h"
#include "evaluation/kuhn_applicator.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::System;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class BS>
class ApplicatorWrapper
  : public ApplicatorInterface<BS>,
    public wrapper< ApplicatorInterface<BS> >
{
  typedef typename BS::real_type R;
  typedef Interval<R> I;
 public:
  ApplicatorWrapper<BS>* clone() const { return this->get_override("clone")(); }
  BS apply(const MapInterface<R>&, const BS&) const {
    return this->get_override("apply")(); }
};


template<class R>
void export_applicator() 
{
  
  class_< ApplicatorWrapper<Rectangle<R> >, boost::noncopyable >("RectangleApplicatorInterface",init<>());
  class_< ApplicatorWrapper<Zonotope<R,ExactTag> >, boost::noncopyable >("ZonotopeApplicatorInterface",init<>());
  class_< ApplicatorWrapper<Zonotope<R,UniformErrorTag> >, boost::noncopyable >("ErrorZonotopeApplicatorInterface",init<>());

  class_< StandardApplicator<R>, 
    bases<ApplicatorInterface< Rectangle<R> >,
          ApplicatorInterface< Zonotope<R,ExactTag> >,
          ApplicatorInterface< Zonotope<R,UniformErrorTag> > > >
    applicator_class("StandardApplicator",init<>());
  applicator_class.def("__call__",(Rectangle<R>(StandardApplicator<R>::*)(const MapInterface<R>&,const Rectangle<R>&)const)&StandardApplicator<R>::apply);
  applicator_class.def("__call__",(Zonotope<R,ExactTag>(StandardApplicator<R>::*)(const MapInterface<R>&,const Zonotope<R,ExactTag>&)const)&StandardApplicator<R>::apply);
  applicator_class.def("__call__",(Zonotope<R,UniformErrorTag>(StandardApplicator<R>::*)(const MapInterface<R>&,const Zonotope<R,UniformErrorTag>&)const)&StandardApplicator<R>::apply);

  class_< KuhnApplicator<R>, bases< ApplicatorInterface< Zonotope<R> > > > kuhn_applicator_class("KuhnApplicator",init<size_type>());
  kuhn_applicator_class.def("__call__",&KuhnApplicator<R>::apply);
  
}

template void export_applicator<FloatPy>();
