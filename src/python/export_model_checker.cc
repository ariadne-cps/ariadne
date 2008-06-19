/***************************************************************************
 *            python/export_model_checker.cc
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

#include "geometry/grid_approximation_scheme.h"
#include "geometry/hybrid_grid_approximation_scheme.h"

#include "system/map.h"
#include "system/vector_field.h"
#include "system/hybrid_automaton.h"

#include "evaluation/evolution_parameters.h"
#include "evaluation/approximator_interface.h"
#include "evaluation/satisfier_interface.h"
#include "evaluation/subdivider_interface.h"
#include "evaluation/applicator_interface.h"
#include "evaluation/integrator_interface.h"
#include "evaluation/evolver_interface.h"
#include "evaluation/map_evolver.h"
#include "evaluation/vector_field_evolver.h"
#include "evaluation/set_based_hybrid_evolver.h"
#include "evaluation/model_checker.h"

using namespace Ariadne;
using namespace Ariadne::Python;

#include <boost/python.hpp>
using namespace boost::python;

template<class T, class Aprx>
void export_system_model_checker()
{
}

template<class R>
void export_model_checker()
{
}

template void export_model_checker<FloatPy>();
