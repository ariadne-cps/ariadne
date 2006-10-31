/***************************************************************************
 *            ariadne.h
 *
 *  Wed Sep 15 15:56 2004
 *  Copyright  2004  Alberto Casagrande, Pieter Collins
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
 
/*! \file ariadne.h
 *  \brief Top-level header file.
 */

#ifndef _ARIADNE_H
#define _ARIADNE_H

#include <gmpxx.h>
#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include <iostream>
#include <iomanip>

/*!
 * \brief Top-level namespace
 */
namespace Ariadne {

/*! \brief Fundamental base types. See the \ref Base module for more information.
 */
namespace Base {}
  
/*! \brief Numerical types, arithmetic, real functions and interval functions. See the \ref Numeric sub-module for more information.
 */
namespace Numeric {}
  
/*! \brief Functions for linear algebra. See the \ref LinearAlgebra sub-module for more information.
 */
namespace LinearAlgebra {}

/*! \brief Combinatorial objects, such as sets and maps on lattices and subdivision trees. See the \ref Combinatoric module for more information.
 */
namespace Combinatoric {}
  
/*! \brief Geometric calculus library. See the \ref Geometry module for more information.
 */
namespace Geometry {}

/*! \brief Classes defining a dynamic system. See the \ref System module for more information.
 */
namespace System {}

/*! \brief Classes for computing system trajectories. See the \ref Evaluation module for more information.
 */
namespace Evaluation {}

/*! \brief Classes for output in various formats. See the \ref Output module for more information.
 */
namespace Output {}

}



#endif /* _ARIADNE_H */
