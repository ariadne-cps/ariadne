/***************************************************************************
 *            geometry.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file geometry.h
 *  \brief Top-level header file for the Geometry module.
 */

#ifndef _ARIADNE_GEOMETRY_H
#define _ARIADNE_GEOMETRY_H

#include <gmpxx.h>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include <iostream>
#include <iomanip>

#include "geometry_declarations.h"

#include "point.h"
#include "rectangle.h"
#include "parallelopiped.h"
#include "polyhedron.h"
#include "simplex.h"

#include "list_set.h"
#include "grid_set.h"
#include "partition_tree_set.h"

#endif /* _ARIADNE_GEOMETRY_H */
