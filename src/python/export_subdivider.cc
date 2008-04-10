/***************************************************************************
 *            python/export_subdivider.cc
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

#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/list_set.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/standard_subdivider.h"

using namespace Ariadne;
using namespace Ariadne::Numeric;
using namespace Ariadne::Geometry;
using namespace Ariadne::Evaluation;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class ES>
class SubdividerWrapper
  : public SubdividerInterface<ES>,
    public wrapper< SubdividerInterface<ES> >
{
  typedef typename ES::real_type R;
  typedef ListSet<ES> ESL;
 public:
  SubdividerWrapper<ES>* clone() const { return this->get_override("clone")(); }
  R radius(const ES&) const { return this->get_override("radius")(); }
  ESL split(const ES&) const { return this->get_override("split")(); }
  ESL subdivide(const ES&, const R&) const { return this->get_override("subdivide")(); }
};


template<class R>
void export_subdivider() 
{
  typedef Zonotope<R> ZES;

  class_< SubdividerWrapper<ZES>, boost::noncopyable >("ZonotopeSubdividerInterface",init<>());

  class_< StandardSubdivider<ZES>, bases< SubdividerInterface<ZES> > >
    subdivider_class("StandardSubdivider",init<>());
}

template void export_subdivider<FloatPy>();
