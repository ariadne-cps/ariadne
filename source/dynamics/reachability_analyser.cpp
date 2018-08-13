/***************************************************************************
 *            reachability_analyser.cpp
 *
 *  Copyright  2006-11  Alberto Casagrande, Pieter Collins, Davide Bresolin
 *
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

#include "../function/functional.hpp"
#include "../config.hpp"

#include <string>
#include <sstream>
#include <algorithm>

#include <list>
#include <set>
#include <vector>
#include <valarray>


#include "../utility/exceptions.hpp"

#include "../numeric/numeric.hpp"

#include "../algebra/vector.hpp"
#include "../algebra/matrix.hpp"

#include "../geometry/box.hpp"
#include "../geometry/list_set.hpp"
#include "../geometry/grid_set.hpp"

#include "../solvers/integrator.hpp"
#include "../solvers/solver.hpp"
#include "../geometry/function_set.hpp"

#include "../dynamics/vector_field.hpp"

#include "../dynamics/orbit.hpp"
#include "../dynamics/vector_field_evolver.hpp"
#include "../dynamics/reachability_analyser.hpp"

#include "../utility/logging.hpp"
#include "../output/graphics.hpp"
#include "../solvers/linear_programming.hpp"

#include "../dynamics/reachability_analyser.tpl.hpp"

namespace Ariadne {

template class ReachabilityAnalyser<VectorField>;
template class ReachabilityAnalyserConfiguration<VectorField>;

} // namespace Ariadne
