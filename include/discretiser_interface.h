/***************************************************************************
 *            discretiser_interface.h
 *
 *  Copyright  2008  Pieter Collins
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

/*! \file discretiser_interface.h
 *  \brief Interface for computing discretised sets and evolution of a system.
 */

#ifndef ARIADNE_DISCRETISER_INTERFACE_H
#define ARIADNE_DISCRETISER_INTERFACE_H

#include "evolver_interface.h"

namespace Ariadne {

using std::pair;


template<class ES> class ListSet;
template<class ES> class Orbit;


//! \ingroup EvaluationModule
//! \brief Interface for evolving a system and discretising the result on a grid.
template<class SYS, class BS>
class DiscretiserInterface
{
  public:
    typedef int AccuracyType;
    typedef SYS SystemType;
    typedef typename SYS::TimeType TimeType;
    typedef BS BasicSetType;

    /*! \brief Destructor. */
    virtual ~DiscretiserInterface() { };

    /*! \brief Make a dynamically-allocated copy. */
    virtual DiscretiserInterface<SystemType,BasicSetType>* clone() const = 0;

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType>
    evolution(const SystemType& system, const BasicSetType& set, const TimeType& time,
              const AccuracyType accuracy, const Semantics semantics) const = 0;

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType>
    lower_evolution(const SystemType& f, const BasicSetType& s, const TimeType& t, const int a) const {
        return evolution(f,s,t,a,LOWER_SEMANTICS); }

    /*! \brief Compute a lower-approximation to the the reachable and evolved sets under the system evolution. */
    virtual Orbit<BasicSetType>
    upper_evolution(const SystemType& f, const BasicSetType& s, const TimeType& t, const int a) const  {
        return evolution(f,s,t,a,UPPER_SEMANTICS); }
};


} // namespace Ariadne



#endif // ARIADNE_EVOLVER_INTERFACE_H
