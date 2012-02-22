/***************************************************************************
 *            paver.h
 *
 *  Copyright  2011-12  Pieter Collins
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

/*! \file paver.h
 *  \brief Class for computing outer approximations to nonlinear sets.
 */

#ifndef ARIADNE_PAVER_H
#define ARIADNE_PAVER_H

#include <iosfwd>

namespace Ariadne {

typedef void Void;
typedef int Int;
typedef std::ostream OutputStream;
class Interval;
template<class X> class Vector;
typedef Vector<Interval> IntervalVector;
template<class X> class VectorFunctionInterface;
typedef VectorFunctionInterface<Interval> IntervalVectorFunctionInterface;
class GridTreeSet;

typedef IntervalVector DomainType;

//! \brief A class for computing outer approximations to sets defined by functions.
class PaverInterface
{
    typedef IntervalVector domain_type;
  public:
    virtual Void
    adjoin_outer_approximation(GridTreeSet& paving, const DomainType& domain, const IntervalVectorFunctionInterface& space_function,
                               const IntervalVectorFunctionInterface& constraint_function, const IntervalVector& constraint_bounds, Int depth) const = 0;

    virtual Void write(OutputStream& os) const = 0;
};

//! \brief A set of the form \f$x=f(s)\f$ for \f$s\in D\f$ satisfying \f$g(s)\leq0\f$ and \f$h(s)=0\f$.
class SubdivisionPaver : public PaverInterface
{
    virtual Void
    adjoin_outer_approximation(GridTreeSet& paving, const DomainType& domain, const IntervalVectorFunctionInterface& space_function,
                               const IntervalVectorFunctionInterface& constraint_function, const IntervalVector& constraint_bounds, Int depth) const;
};

//! \brief A set of the form \f$x=f(s)\f$ for \f$s\in D\f$ satisfying \f$g(s)\leq0\f$ and \f$h(s)=0\f$.
class OptimalConstraintPaver : public PaverInterface
{
    virtual Void
    adjoin_outer_approximation(GridTreeSet& paving, const DomainType& domain, const IntervalVectorFunctionInterface& space_function,
                               const IntervalVectorFunctionInterface& constraint_function, const IntervalVector& constraint_bounds, Int depth) const;
};

} //namespace Ariadne

#endif /* ARIADNE_PAVER_H */
