/***************************************************************************
 *            numerical_types.h
 *
 *  Copyright  2004-6  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
/*! \file numerical_types.h
 *  \brief Type definitions and conversion operators for fundamental %Ariadne types.
 */

#ifndef _ARIADNE_NUMERICAL_TYPES_H
#define _ARIADNE_NUMERICAL_TYPES_H

#ifdef DOUBLE_REAL 
#include "float64.h"
#endif

#ifdef MPFLOAT_REAL
#include "mpfloat.h"
#endif

#ifdef DYADIC_REAL
#include "dyadic.h"
#endif

#ifdef RATIONAL_REAL
#include "rational.h"
#endif

#endif /* _ARIADNE_NUMERICAL_TYPES */
