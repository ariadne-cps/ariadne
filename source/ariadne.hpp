/***************************************************************************
 *            ariadne.hpp
 *
 *  Copyright 2008  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file ariadne.hpp
 *  \brief Top-level header file includes all user headers.
 */

#ifndef ARIADNE_ARIADNE_HPP
#define ARIADNE_ARIADNE_HPP

//! \brief Top-level %Ariadne namespace
namespace Ariadne {
}

#include <cstdlib>
#include <iostream>
#include <iomanip>

#include "config.hpp"

#include "numeric/numeric.hpp"
#include "numeric/rational.hpp"
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"

#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "function/taylor_function.hpp"

#include "set/function_set.hpp"
#include "set/affine_set.hpp"
#include "paving/grid_paving.hpp"

#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "set/curve.hpp"
#include "dynamics/orbit.hpp"

#include "integration/integrator.hpp"
#include "equation-solving/equation_solver.hpp"

#include "hybrid/discrete_location.hpp"
#include "hybrid/discrete_event.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "evolution/hybrid_time.hpp"
#include "hybrid/hybrid_automata.hpp"

#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/differential_inclusion_evolver.hpp"
#include "dynamics/vector_field_simulator.hpp"
#include "dynamics/iterated_map_evolver.hpp"
#include "evolution/hybrid_evolver.hpp"
#include "evolution/hybrid_simulator.hpp"
#include "hybrid/hybrid_reachability_analyser.hpp"

#include "labelled/figure.hpp"
#include "rendering/rendering_manager.hpp"
#include "drawing/drawer_manager.hpp"
// Cairo-specific interfaces are no longer part of the extracted rendering API.
#include "drawing/drawing_command_line_interface.hpp"
#include "hybrid/hybrid_graphics.hpp"

#include "logging/logging.hpp"
#include "utility/stlio.hpp"

#endif
