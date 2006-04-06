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

typedef RRectangle (*MapRectBinFun) (const RMapBase&, const RRectangle&);
typedef RParallelotope (*MapPltpBinFun) (const RMapBase&, const RParallelotope&);
typedef RParallelotopeListSet (*MapLSPltpBinFun) (const RMapBase&, const RParallelotopeListSet&);
typedef RGridMaskSet (*CRFun) (const RMapBase&, const RRectangleListSet&, const RFiniteGrid&, const RRectangle&);

void export_apply() {
  def("apply", MapRectBinFun(&apply), "apply the image of a map to a set" );
  def("apply", MapPltpBinFun(&apply), "apply the image of a map to a set" );
  def("apply", MapLSPltpBinFun(&apply), "apply the image of a map to a set" );
  def("chainreach", CRFun(&chainreach), "chain reach of a set" );
}
