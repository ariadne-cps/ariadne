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

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

#include "real.h"
#include "function.h"
#include "user_function.h"
#include "taylor_model.h"

#include "space.h"
#include "expression.h"
#include "variables.h"
#include "assignment.h"

#include "function_set.h"
#include "taylor_set.h"
#include "grid_set.h"
#include "hybrid_set.h"

#include "point.h"
#include "box.h"
#include "polytope.h"
#include "polyhedron.h"

#include "hybrid_time.h"
#include "orbit.h"

#include "discrete_state.h"
#include "discrete_event.h"
#include "hybrid_automaton.h"
#include "hybrid_io_automaton.h"

#include "vector_field_evolver.h"
#include "hybrid_evolver.h"
#include "discretiser.h"
#include "parametric.h"
#include "reachability_analyser.h"
#include "verifier.h"

#include "serialization.h"
#include "graphics.h"
#include "textplot.h"

#endif
