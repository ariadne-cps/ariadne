/***************************************************************************
 *            python/export_reducer.cc
 *
 *  Copyright  2008  Pieter Collins
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

#include "geometry/box.h"
#include "geometry/zonotope.h"
#include "geometry/taylor_set.h"
#include "evaluation/reducer_interface.h"
#include "evaluation/identity_reducer.h"
#include "evaluation/orthogonal_reducer.h"
#include "evaluation/cascade_reducer.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class ES>
class ReducerWrapper
  : public ReducerInterface<ES>,
    public wrapper< ReducerInterface<ES> >
{
 public:
  ReducerWrapper<ES>* clone() const { return this->get_override("clone")(); }
  ES over_approximate(const ES&) const { return this->get_override("over_approximate")(); }
};


template<class R>
void export_reducer() 
{
  typedef Zonotope<R> ZES;

  class_< ReducerWrapper<ZES>, boost::noncopyable >("ZonotopeReducerInterface",init<>());

  class_< IdentityReducer<ZES>, bases< ReducerInterface<ZES> > >
    identity_reducer_class("IdentityReducer",init<>());
  class_< CascadeReducer<ZES>, bases< ReducerInterface<ZES> > >
    cascade_reducer_class("CascadeReducer",init<uint>());
  class_< OrthogonalReducer<ZES>, bases< ReducerInterface<ZES> > >
    orthogonal_reducer_class("OrthogonalReducer",init<>());
}

template void export_reducer<FloatPy>();
