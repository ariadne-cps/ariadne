/***************************************************************************
 *            numeric/float.h
 *
 *  Copyright  2007 Pieter Collins
 *  pieter.collins@cwi.nl
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
 
/*! \file numeric/float.h
 *  \brief Type definitions and conversion operators for floating point numbers.
 */

#ifndef ARIADNE_NUMERIC_FLOAT_H
#define ARIADNE_NUMERIC_FLOAT_H

#include "config.h"

#include "numeric/float.template.h"

#ifdef ENABLE_FLOAT64
  #include "numeric/float64.h"
#endif

#ifdef ENABLE_FLOATMP
  #include "numeric/floatmp.h"
#endif

#include "numeric/arithmetic.h"

#endif /* ARIADNE_NUMERIC_FLOAT_H */
