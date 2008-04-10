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
#include "system/map.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/standard_applicator.h"

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
  BS apply(const Map<R>&, const BS&) const {
    return this->get_override("apply")(); }
};


template<class R>
class Applicator 
  : public StandardApplicator< Rectangle<R> >,
    public StandardApplicator< Zonotope<R> >
{ };



template<class R>
void export_applicator() 
{
  
  class_< ApplicatorWrapper<Rectangle<R> >, boost::noncopyable >("RectangleApplicatorInterface",init<>());
  class_< ApplicatorWrapper<Zonotope<R> >, boost::noncopyable >("ZonotopeApplicatorInterface",init<>());

  class_< Applicator<R>, 
    bases<ApplicatorInterface< Rectangle<R> >,
          ApplicatorInterface< Zonotope<R> > > >
    applicator_class("StandardApplicator",init<>());
  applicator_class.def("__call__",(Rectangle<R>(Applicator<R>::*)(const Map<R>&,const Rectangle<R>&)const)&StandardApplicator< Rectangle<R> >::apply);
  applicator_class.def("__call__",(Zonotope<R>(Applicator<R>::*)(const Map<R>&,const Zonotope<R>&)const)&StandardApplicator< Zonotope<R> >::apply);

  
}

template void export_applicator<FloatPy>();
