/***************************************************************************
 *            sweeper.h
 *
 *  Copyright 2010-11  Pieter Collins
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

/*! \file sweeper.h
 *  \brief Classes for simplifying the representation of a polynomial expansion.
 */

#ifndef ARIADNE_SWEEPER_H
#define ARIADNE_SWEEPER_H

#include "utility/macros.h"
#include "utility/attribute.h"
#include "algebra/multi_index.h"
#include "algebra/expansion.h"

namespace Ariadne {

class Float;
template<class X> class Expansion;

template<class SWP> class SweeperBase;

class SweeperInterface {
    friend class Sweeper;
  public:
    inline bool discard(const MultiIndex& a, const Float& x) const { return this->_discard(a,x); }
    inline void sweep(Expansion<Float>& p, Float& e) const { this->_sweep(p,e); }
    inline void sweep(Expansion<Float>& p) const { this->_sweep(p); }
  private:
    virtual SweeperInterface* _clone() const = 0;
    virtual void _sweep(Expansion<Float>& p, Float& e) const = 0;
    virtual void _sweep(Expansion<Float>& p) const = 0;
    virtual bool _discard(const MultiIndex& a, const Float& x) const = 0;
    virtual void _write(std::ostream& os) const = 0;
    friend std::ostream& operator<<(std::ostream& os, const SweeperInterface& swp) { swp._write(os); return os; }
};

/*! \brief A class for removing terms in a polynomial expansion which fall below some significance level, putting them in a uniform error.
 *
 *  \internal The virtual methods in SweeperInterface are private so that they can be implemented in terms of public inline methods for efficiency.
 *  In particular, the concrete sweeper classes have an inline discard(...) method which is used in the sweep(...) method since the exact type is known.
 *  This is a major efficiency improvement.
 */
class Sweeper {
  public:
    Sweeper();
    inline Sweeper(const SweeperInterface& p) : _ptr(p._clone()) { }
    inline Sweeper(std::shared_ptr<const SweeperInterface> p) : _ptr(p) { }
    inline Sweeper(const SweeperInterface* p) : _ptr(p) { }
    inline operator const SweeperInterface& () const { return *_ptr; }
    inline std::shared_ptr<const SweeperInterface> operator& () { return _ptr; }
  public:
    //! \brief Returns \a true if the term with index \a a and coefficient \a x should be discarded.
    inline bool discard(const MultiIndex& a, const Float& x) const { return this->_ptr->_discard(a,x); }
    //! \brief Discard terms in the expansion, adding the absolute value of the coefficient to the uniform error.
    inline void sweep(Expansion<Float>& p, Float& e) const { this->_ptr->_sweep(p,e); }
    //! \brief Discard terms in the expansion, without keeping track of discarded terms.
    inline void sweep(Expansion<Float>& p) const { this->_ptr->_sweep(p); }
    friend std::ostream& operator<<(std::ostream& os, const Sweeper& swp) { return os << *swp._ptr; }
  private:
    std::shared_ptr<const SweeperInterface> _ptr;
};

template<class SWP> class SweeperBase
    : public virtual SweeperInterface
{
    virtual SweeperInterface* _clone() const;
    virtual bool _discard(const MultiIndex& a, const Float& x) const;
    virtual void _sweep(Expansion<Float>& p, Float& e) const;
    virtual void _sweep(Expansion<Float>& p) const;
};


template<class SWP>
SweeperInterface* SweeperBase<SWP>::_clone() const
{
    return new SWP(static_cast<const SWP&>(*this));
}

template<class SWP>
bool SweeperBase<SWP>::_discard(const MultiIndex& a, const Float& x) const
{
    return static_cast<const SWP*>(this)->SWP::discard(a,x);
}

template<class SWP>
void SweeperBase<SWP>::_sweep(Expansion<Float>& p, Float& e) const
{
    Expansion<Float>::const_iterator end=p.end();
    Expansion<Float>::const_iterator adv=p.begin();
    Expansion<Float>::iterator curr=p.begin();
    set_rounding_upward();
    Float te=0.0;
    while(adv!=end) {
        if(static_cast<const SWP*>(this)->SWP::discard(adv->key(),adv->data())) {
            te+=abs(adv->data());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;
    p.resize(curr-p.begin());
    set_rounding_to_nearest();
}

template<class SWP>
void SweeperBase<SWP>::_sweep(Expansion<Float>& p) const
{
    Expansion<Float>::const_iterator end=p.end();
    Expansion<Float>::const_iterator adv=p.begin();
    Expansion<Float>::iterator curr=p.begin();
    while(adv!=end) {
        if(static_cast<const SWP*>(this)->SWP::discard(adv->key(),adv->data())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(curr-p.begin());
}

//! \brief A sweeper class which discards terms whose absolute value is smaller than a threshold.
class ThresholdSweeper : public SweeperBase<ThresholdSweeper> {
    double _sweep_threshold;
  public:
    ThresholdSweeper(double sweep_threshold) : _sweep_threshold(sweep_threshold) { ARIADNE_ASSERT(sweep_threshold>=0.0); }
    Float sweep_threshold() const { return _sweep_threshold; }
    inline bool discard(const MultiIndex& a, const Float& x) const { return abs(x) < this->_sweep_threshold; }
  private:
    virtual void _write(std::ostream& os) const { os << "ThresholdSweeper( sweep_threshold="<<this->_sweep_threshold<<" )"; };
};

//! \brief A sweeper class which does not discard any terms at all.
class TrivialSweeper : public SweeperBase<TrivialSweeper> {
  public:
    inline bool discard(const MultiIndex& a, const Float& x) const { return false; }
  private:
    virtual void _sweep(Expansion<Float>& p, Float& e) const { }
    virtual void _write(std::ostream& os) const { os << "TrivialSweeper"; }
};

//! \brief A sweeper class which only discards the zero term.
class NullSweeper : public SweeperBase<NullSweeper> {
  public:
    inline bool discard(const MultiIndex& a, const Float& x) const { return x==0.0; }
  private:
    virtual void _write(std::ostream& os) const { os << "NullSweeper"; }
};

//! \brief A sweeper class which discards non-affine terms.
class AffineSweeper : public SweeperBase<AffineSweeper> {
  public:
    inline bool discard(const MultiIndex& a, const Float& x) const { return a.degree()>1; }
  private:
    virtual void _write(std::ostream& os) const { os << "AffineSweeper"; }
};

//! \brief A sweeper class which discards terms whose total degree is above some threshold.
class GradedSweeper : public SweeperBase<GradedSweeper> {
  public:
    GradedSweeper(uint degree) : _degree(degree) { }
    uint degree() const { return this->_degree; }
    inline bool discard(const MultiIndex& a, const Float& x) const { return a.degree()>this->_degree; }
  private:
    virtual void _write(std::ostream& os) const { os << "GradedSweeper( degree="<<this->_degree<<" )"; }
  private:
    uint _degree;
};

inline Sweeper::Sweeper() : _ptr(new ThresholdSweeper(std::numeric_limits<float>::epsilon())) { }


} // namespace Ariadne

#endif // ARIADNE_SWEEPER_H
