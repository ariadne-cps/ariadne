/***************************************************************************
 *            hybrid_time.h
 *
 *  Copyright 2008-17  Pieter Collins
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

#include "numeric/numeric.h"

namespace Ariadne {

typedef Integer DiscreteTimeType;

//! \ingroup SystemModule
//! \brief A value in a hybrid time domain, being a pair comprising a real \a continuous_time
//! and an integer \a discrete_time.
//!
//! When a %HybridTime is used to define a bound on a hybrid evolution, the evolution should
//! stop when <em>either</em> the continuous time or the discrete time reaches the bounding
//! value. This is to ensure that the evolution time is finite; in particular, that no
//! Zeno behaviour occurs.
class HybridTime
{
  public:
    //! \brief The type used for continuous (real, physical) time.
    typedef Real ContinuousTimeType;
    //! \brief The type used for discrete time (steps).
    typedef Integer DiscreteTimeType;

    //! \brief The continuous (real, physical) time.
    const ContinuousTimeType& continuous_time() const { return this->_continuous_time; }
    //! \brief The number of discrete steps taken.
    const DiscreteTimeType& discrete_time() const { return this->_discrete_time; }
  private: public:
    //! \brief The continuous (real, physical) time.
    ContinuousTimeType _continuous_time;
    //! \brief The number of discrete steps taken.
    DiscreteTimeType _discrete_time;
  public:
    HybridTime(Real t, Integer n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(RawFloat64 t, Integer n)
      : _continuous_time(Float64Value(t)), _discrete_time(n) { }
    HybridTime(Float64Value t, Integer n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(double t, Int n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(Int n, double t)
      : _continuous_time(t), _discrete_time(n) {
          ARIADNE_FAIL_MSG("HybridTime(Int,double) is incorrect; use HybridTime(Real,Integer) instead."); }
};

inline Bool same(const HybridTime& ht1, const HybridTime& ht2) {
    return same(ht1._continuous_time,ht2._continuous_time) &&
        ht1._discrete_time==ht2._discrete_time;
}

inline ValidatedNegatedSierpinskian operator==(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1._continuous_time==ht2._continuous_time &&
        ht1._discrete_time==ht2._discrete_time;
}

inline ValidatedSierpinskian operator!=(const HybridTime& ht1, const HybridTime& ht2) {
    return ht1._continuous_time!=ht2._continuous_time ||
        ht1._discrete_time!=ht2._discrete_time;
}

inline ValidatedKleenean operator<=(const HybridTime& ht1, const HybridTime& ht2) {
    return ValidatedKleenean(ht1._continuous_time<=ht2._continuous_time) &&
        Boolean(ht1._discrete_time<=ht2._discrete_time);
}

inline ValidatedKleenean operator<(const HybridTime& ht1, const HybridTime& ht2) {
    return ValidatedKleenean(ht1._continuous_time< ht2._continuous_time) &&
        Boolean(ht1._discrete_time<=ht2._discrete_time);
}

inline OutputStream& operator<<(OutputStream& os, const HybridTime& ht) {
    return os << "("<<ht._continuous_time<<","<<ht._discrete_time<<")";
}

} // namespace Ariadne

#endif // ARIADNE_HYBRID_TIME_H
