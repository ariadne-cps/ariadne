/***************************************************************************
 *            hybrid_time.h
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

/*! \file hybrid_time.h
 *  \brief Hybrid times
 */

#ifndef ARIADNE_HYBRID_TIME_H
#define ARIADNE_HYBRID_TIME_H

#include "numeric.h"

namespace Ariadne {

//! \ingroup SystemModule
//! \brief A value in a hybrid time domain, being a pair comprising a real \a continuous_time
//! and an integer \a discrete_time.
//!
//! When a %HybridTime is used to define a bound on a hybrid evolution, the evolution should
//! stop when <em>either</em> the continuous time or the discrete time reaches the bounding
//! value. This is to ensure that the evolution time is finite; in particular, that no
//! Zeno behaviour occurs.
struct HybridTime
{
    //! \brief The type used for continuous (real, physical) time.
    typedef Real ContinuousTimeType;
    //! \brief The type used for discrete time (steps).
    typedef int DiscreteTimeType;

    //! \brief The continuous (real, physical) time.
    const ContinuousTimeType& continuous_time() const { return this->_continuous_time; }
    //! \brief The number of discrete steps taken.
    const DiscreteTimeType& discrete_time() const { return this->_discrete_time; }

    //! \brief The continuous (real, physical) time.
    ContinuousTimeType _continuous_time;
    //! \brief The number of discrete steps taken.
    DiscreteTimeType _discrete_time;
  public:
    HybridTime(Real t, int n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(Float t, int n)
      : _continuous_time(Dyadic(t)), _discrete_time(n) { }
    HybridTime(Dyadic t, int n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(double t, int n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(int n, double t)
      : _continuous_time(t), _discrete_time(n) {
          ARIADNE_FAIL_MSG("HybridTime(int,double) is incorrect; use HybridTime(Real,Integer) instead."); }
};

inline bool operator==(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1._continuous_time==ht2._continuous_time &&
        ht1._discrete_time==ht2._discrete_time;
}

inline bool operator!=(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1._continuous_time!=ht2._continuous_time ||
        ht1._discrete_time!=ht2._discrete_time;
}

inline bool operator<=(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1._continuous_time<=ht2._continuous_time &&
        ht1._discrete_time<=ht2._discrete_time;
}

inline bool operator<(const HybridTime& ht1, const HybridTime& ht2) {
    return (ht1<=ht2) && (ht1 != ht2);
}

inline std::ostream& operator<<(std::ostream& os, const HybridTime& ht) {
    return os << "("<<ht._continuous_time<<","<<ht._discrete_time<<")";
}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_TIME_H
