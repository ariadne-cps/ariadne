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

/*! \class interval
 * \brief A templated class representing an interval of real numbers.
 *
 * An interval of real numbers with endpoints of type \a R.
 * All operations on an interval must be guarenteed to return an interval contining the exact result.
 * If \a val of real numbers with endpoints of type \a R.
 * All operations on an interval must be guarenteed to return an interval contining the exact result.
 * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
 * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a maximum error of computation.
 *
 * COMMENT FOR ALBERTO: How should we specify error bounds for computations on types which support arbitrary-precision computing?
 * I would suggest using a global(ish) "precision" variable which either may be "locked" by a computation, or stored on a stack.
 * I think functions should keep a natural syntax, so we can declare <code>rational cos(rational)</code>,
 * and set error bounds for the computation elsewhere.
 * Of course, this partly depends on the interval / rational arithmetic library we use.
 *
 * Currently implemented using the boost::numeric::interval from the Boost C++ library.
 */
using boost::numeric::interval;



/*! \brief Geometric calculus library.
 */
namespace Geometry {}

/*! \brief Classes defining a hybrid system.
 */
namespace HybridDefinitions {}

/*! \brief Classes defining map.
 */
namespace Map {

/*! \brief Classes defining affine map.
 */
namespace Affine {}

}

/*! \brief Classes defining vector field.
 *
 *  In this namespace there are also vector field integrator.
 */
namespace VectorField {

/*! \brief Classes defining affine vector field.
 *
 * In this namespace there are also affine vector field integrator.
 */
namespace Affine {}

}

/*! \brief Functions for linear algebra.
 */
namespace LinearAlgebra {}


/*! \brief Functions for computing hybrid system trajectories.
 */
namespace Evaluation {}

}

/* No input routine for intervals defined by boost */
namespace boost {
    namespace numeric {
	template<typename R> 
	std::istream&
	operator>>(std::istream& is, interval<R>& ivl) 
	{
	    R l,u;
	    char c1,c2,c3;
	    is >> c1 >> l >> c2 >> u >> c3;
	    ivl=interval<R>(l,u);
	}
    }
}

#endif /* _ARIADNE_H */
