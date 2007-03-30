/***************************************************************************
 *            python/python_float.h
 *
 *  Copyright  2005  Alberto Casagrande, Pieter Collins
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

/*! \file python_utilities.h
 *  Commonly used inline methods for the Python interface.
 */
 
#ifndef ARIADNE_PYTHON_FLOAT_H
#define ARIADNE_PYTHON_FLOAT_H

#include <config.h>

#if PYTHON_FLOAT == Float64 
#include "numeric/float64.h" 
namespace Ariadne { typedef Numeric::Float64 Float; }
#elif PYTHON_FLOAT == FloatMP 
#include "numeric/floatmp.h" 
namespace Ariadne { typedef Numeric::FloatMP Float; }
#endif

#endif /* ARIADNE_PYTHON_FLOAT_H */
