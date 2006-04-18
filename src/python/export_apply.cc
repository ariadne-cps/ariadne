/***************************************************************************
 *            python/export_apply.cc
 *
 *  6 February 2006
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
 ****************************************************************************/

/*
 *  This program is free software; you can rediself_ns::stribute it and/or modify
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

#include "geometry/rectangle.h"
#include "geometry/parallelotope.h"
#include "geometry/list_set.h"
#include "geometry/grid_set.h"
#include "evaluation/apply.h"

#include "python/typedefs.h"
using namespace Ariadne;
using namespace Ariadne::Evaluation;

#include <boost/python.hpp>
using namespace boost::python;

typedef C0Applicator<Real> RC0Applicator;
typedef C1Applicator<Real> RC1Applicator;

typedef RRectangle (RC0Applicator::*ApplMapRectBinFunc) (const RMapBase&, const RRectangle&) const;
typedef RParallelotope (RC1Applicator::*ApplMapPltpBinFunc) (const RMapBase&, const RParallelotope&) const;
typedef RParallelotopeListSet (RC1Applicator::*ApplMapLSPltpBinFunc) (const RMapBase&, const RParallelotopeListSet&) const;
typedef RGridMaskSet (RC1Applicator::*ApplMapGMSFunc) (const RMapBase&, const RGridMaskSet&, const RGridMaskSet&) const;

typedef RRectangle (*MapRectBinFunc) (const RMapBase&, const RRectangle&);
typedef RParallelotope (*MapPltpBinFunc) (const RMapBase&, const RParallelotope&);
typedef RParallelotopeListSet (*MapLSPltpBinFunc) (const RMapBase&, const RParallelotopeListSet&);
typedef RGridMaskSet (*MapGMSFunc) (const RMapBase&, const RGridMaskSet&, const RGridMaskSet&);

void export_apply() {

  class_<RC1Applicator>("C1Applicator",init<>())
    .def("apply", ApplMapRectBinFunc(&RC0Applicator::apply), "apply the image of a map to a set" )
    .def("apply", ApplMapPltpBinFunc(&RC1Applicator::apply), "apply the image of a map to a set" )
    .def("apply", ApplMapLSPltpBinFunc(&RC1Applicator::apply), "apply the image of a map to a set" )
    .def("apply", ApplMapGMSFunc(&RC1Applicator::apply), "apply the image of a map to a set" )
    .def("chainreach", ApplMapGMSFunc(&RC1Applicator::chainreach), "Compute the chain reachable set")
  ;
  
  def("apply", MapRectBinFunc(&apply), "apply the image of a map to a set" );
  def("apply", MapPltpBinFunc(&apply), "apply the image of a map to a set" );
  def("apply", MapLSPltpBinFunc(&apply), "apply the image of a map to a set" );
  def("apply", MapGMSFunc(&apply), "apply the image of a map to a set" );
  def("chainreach", MapGMSFunc(&chainreach), "Compute the chain reachable set");
  
}
