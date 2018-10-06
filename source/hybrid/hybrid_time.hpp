/***************************************************************************
 *            hybrid_time.hpp
 *
 *  Copyright 2008-17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

/*! \file hybrid_time.hpp
 *  \brief Hybrid times
 */

#ifndef ARIADNE_HYBRID_TIME_HPP
#define ARIADNE_HYBRID_TIME_HPP

#include "../numeric/numeric.hpp"

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
    explicit HybridTime(Real t)
      : _continuous_time(t), _discrete_time(0) { }
    HybridTime(Real t, Integer n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(RawFloatDP t, Integer n)
      : _continuous_time(FloatDPValue(t)), _discrete_time(n) { }
    HybridTime(FloatDPValue t, Integer n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(double t, Int n)
      : _continuous_time(t), _discrete_time(n) { }
    HybridTime(Int n, double t)
      : _continuous_time(t), _discrete_time(n) {
          ARIADNE_FAIL_MSG("HybridTime(Int,double) is incorrect; use HybridTime(Real,Integer) instead."); }

    friend HybridTime operator*(const Integer& c, const HybridTime& ht) {
        return HybridTime(c*ht._continuous_time,c*ht._discrete_time);
    }

    friend HybridTime operator-(const HybridTime& ht1, const HybridTime& ht2) {
        return HybridTime(ht1._continuous_time-ht2._continuous_time,max(ht1._discrete_time-ht2._discrete_time,Integer(0)));
    }

    friend Bool same(const HybridTime& ht1, const HybridTime& ht2) {
        return same(ht1._continuous_time,ht2._continuous_time) &&
            ht1._discrete_time==ht2._discrete_time;
    }

    friend NegatedSierpinskian operator==(const HybridTime& ht1, const HybridTime& ht2) {
        return ht1._continuous_time==ht2._continuous_time &&
            ht1._discrete_time==ht2._discrete_time;
    }

    friend Sierpinskian operator!=(const HybridTime& ht1, const HybridTime& ht2) {
        return ht1._continuous_time!=ht2._continuous_time ||
            ht1._discrete_time!=ht2._discrete_time;
    }

    friend Kleenean operator<=(const HybridTime& ht1, const HybridTime& ht2) {
        return Kleenean(ht1._continuous_time<=ht2._continuous_time) &&
            Boolean(ht1._discrete_time<=ht2._discrete_time);
    }

    friend Kleenean operator<(const HybridTime& ht1, const HybridTime& ht2) {
        return Kleenean(ht1._continuous_time< ht2._continuous_time) &&
            Boolean(ht1._discrete_time<=ht2._discrete_time);
    }

    friend Kleenean operator>(const HybridTime& ht1, const HybridTime& ht2) {
        return Kleenean(ht1._continuous_time> ht2._continuous_time) &&
            Boolean(ht1._discrete_time>=ht2._discrete_time);
    }

    friend Kleenean operator>(const HybridTime& ht1, const ContinuousTimeType& ct2) {
        return Kleenean(ht1._continuous_time> ct2);
    }

    friend OutputStream& operator<<(OutputStream& os, const HybridTime& ht) {
        return os << "("<<ht._continuous_time<<","<<ht._discrete_time<<")";
    }
};


} // namespace Ariadne

#endif // ARIADNE_HYBRID_TIME_HPP
