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
//     typedef boost::numeric::interval interval;

/*! \brief An integer
 * 
 * An element of the ring of integers.
 * Must allow denotation of any integer, including arbitrarily large values.
 * Integer quotient and remainder must be supported.
 *
 * Currently implemented using mpz_class from the GNU Multiple Precision Library.
 */
typedef mpz_class integer;
 
/*! A dyadic rational (i.e. of form @a m/2^n).
 * 
 * An element of the ring of dyadic rationals.
 * Must allow denotation of any dyadic rational.
 * May be created without loss of precision from any integral or floating point type, 
 * or from any rational of the form m/2^n.
 *
 * Currently implemented using mpf_class from the GNU Multiple Precision library.
 */
typedef mpf_class dyadic;
 
/*! \brief A rational number.
 * 
 * An element of the field of rationals.
 * Must allow denotation of any rational.
 * May be created without loss of precision from any integral or floating point type, and from a dyadic.
 *
 * Currently implemented using mpq_class from the GNU Multiple Precision library.
 */
typedef mpq_class rational;
 

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


namespace boost {
namespace numeric {
    template<class Ch, class ChTr>
    std::basic_ostream<Ch, ChTr>&
    operator<<(std::basic_ostream<Ch, ChTr>& os, const interval<Ariadne::rational>& r)
    {
	typename std::basic_ostream<Ch, ChTr>::sentry sentry(os);
	if (sentry) {
	    Ariadne::rational l = r.lower(), u = r.upper();
	    os << '[' << l << ',' << u << ']';
	}
	return os;
    }
}
}

std::ostream& operator<<(std::ostream& os, const boost::numeric::interval<Ariadne::rational>& r) {
    os << "[" << r.lower() << "," << r.upper() << "]";
    return os;
}

#endif /* _ARIADNE_H */
