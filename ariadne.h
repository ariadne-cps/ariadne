/***************************************************************************
 *            ariadne.h
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

/*! \file ariadne.h
 *  \brief Top-level header file includes all user headers.
 */

#ifndef ARIADNE_ARIADNE_H
#define ARIADNE_ARIADNE_H

//! \brief Top-level %Ariadne namespace
namespace Ariadne {
}

#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>

using std::cin; using std::cout; using std::cerr; using std::clog;
using std::endl; using std::flush;
using std::ofstream; using std::ifstream;

#include "config.h"

#include "numeric/numeric.h"
#include "algebra/vector.h"
#include "algebra/matrix.h"

#include "function/function.h"
#include "function/taylor_model.h"

#include "geometry/function_set.h"
#include "geometry/affine_set.h"
#include "geometry/grid_set.h"

#include "geometry/point.h"
#include "geometry/box.h"
#include "geometry/curve.h"
#include "dynamics/orbit.h"

#include "solvers/integrator.h"
#include "solvers/solver.h"

#include "hybrid/discrete_location.h"
#include "hybrid/discrete_event.h"
#include "hybrid/hybrid_set.h"
#include "hybrid/hybrid_orbit.h"
#include "hybrid/hybrid_time.h"
#include "hybrid/hybrid_automaton.h"

#include "solvers/vector_field_evolver.h"
#include "hybrid/hybrid_evolver.h"
#include "hybrid/hybrid_reachability_analyser.h"

#ifdef ARIADNE_ENABLE_SERIALIZATION
#include "output/serialization.h"
#endif /* ARIADNE_ENABLE_SERIALIZATION */

#include "output/graphics.h"
#include "hybrid/hybrid_graphics.h"

#endif
