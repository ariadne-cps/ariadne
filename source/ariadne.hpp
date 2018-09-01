/***************************************************************************
 *            ariadne.hpp
 *
 *  Copyright 2008  Pieter Collins
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
#include "algebra/vector.hpp"
#include "algebra/matrix.hpp"

#include "function/function.hpp"
#include "function/taylor_model.hpp"

#include "geometry/function_set.hpp"
#include "geometry/affine_set.hpp"
#include "geometry/grid_paving.hpp"

#include "geometry/point.hpp"
#include "geometry/box.hpp"
#include "geometry/curve.hpp"
#include "dynamics/orbit.hpp"

#include "solvers/integrator.hpp"
#include "solvers/solver.hpp"

#include "hybrid/discrete_location.hpp"
#include "hybrid/discrete_event.hpp"
#include "hybrid/hybrid_set.hpp"
#include "hybrid/hybrid_paving.hpp"
#include "hybrid/hybrid_orbit.hpp"
#include "hybrid/hybrid_time.hpp"
#include "hybrid/hybrid_automata.hpp"

#include "dynamics/vector_field_evolver.hpp"
#include "dynamics/differential_inclusion.hpp"
#include "hybrid/hybrid_evolver.hpp"
#include "hybrid/hybrid_simulator.hpp"
#include "hybrid/hybrid_reachability_analyser.hpp"

#include "output/graphics.hpp"
#include "hybrid/hybrid_graphics.hpp"

#endif
