/***************************************************************************
 *            sweeper.hpp
 *
 *  Copyright 2010-17  Pieter Collins
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

/*! \file sweeper.hpp
 *  \brief Classes for simplifying the representation of a polynomial expansion.
 */

#ifndef ARIADNE_SWEEPER_HPP
#define ARIADNE_SWEEPER_HPP

#include "utility/macros.hpp"
#include "utility/attribute.hpp"
#include "numeric/float.decl.hpp"
#include "algebra/multi_index.hpp"
#include "algebra/expansion.hpp"

namespace Ariadne {

template<class F> class Sweeper;
template<class X> class Expansion;

template<class F> class SweeperInterface {
    friend class Sweeper<F>;
  protected:
    typedef typename F::PrecisionType PR;
  public:
    inline Bool discard(const MultiIndex& a, const F& x) const { return this->_discard(a,x); }
    inline Void sweep(Expansion<F>& p, F& e) const { this->_sweep(p,e); }
    inline Void sweep(Expansion<F>& p) const { this->_sweep(p); }
    inline PR precision() const { return this->_precision(); }
  private:
    virtual SweeperInterface* _clone() const = 0;
    virtual PR _precision() const = 0;
    virtual Void _sweep(Expansion<F>& p, F& e) const = 0;
    virtual Void _sweep(Expansion<F>& p) const = 0;
    virtual Bool _discard(const MultiIndex& a, const F& x) const = 0;
    virtual Void _write(OutputStream& os) const = 0;
    friend OutputStream& operator<<(OutputStream& os, const SweeperInterface& swp) { swp._write(os); return os; }
};

/*! \brief A class for removing terms in a polynomial expansion which fall below some significance level, putting them in a uniform error.
 *
 *  \internal The virtual methods in SweeperInterface are private so that they can be implemented in terms of public inline methods for efficiency.
 *  In particular, the concrete sweeper classes have an inline discard(...) method which is used in the sweep(...) method since the exact type is known.
 *  This is a major efficiency improvement.
 */
template<class F> class Sweeper {
    typedef typename F::PrecisionType PR;
  public:
    typedef F FloatType;
    typedef PR PrecisionType;
    typedef MultiIndex IndexType;
  public:
    Sweeper();
    inline Sweeper(const SweeperInterface<F>& p) : _ptr(p._clone()) { }
    inline Sweeper(std::shared_ptr<const SweeperInterface<F>> p) : _ptr(p) { }
    inline Sweeper(const SweeperInterface<F>* p) : _ptr(p) { }
    inline operator const SweeperInterface<F>& () const { return *_ptr; }
    inline std::shared_ptr<const SweeperInterface<F>> operator& () { return _ptr; }
  public:
    //! \brief The precision to which terms should be built.
    inline PrecisionType precision() const { return this->_ptr->_precision(); }
    //! \brief Returns \a true if the term with index \a a and coefficient \a x should be discarded.
    inline Bool discard(const IndexType& a, const FloatType& x) const { return this->_ptr->_discard(a,x); }
    //! \brief Discard terms in the expansion, adding the absolute value of the coefficient to the uniform error.
    inline Void sweep(Expansion<FloatType>& p, FloatType& e) const { this->_ptr->_sweep(p,e); }
    //! \brief Discard terms in the expansion, without keeping track of discarded terms.
    inline Void sweep(Expansion<FloatType>& p) const { this->_ptr->_sweep(p); }
    friend OutputStream& operator<<(OutputStream& os, const Sweeper<F>& swp) { return os << *swp._ptr; }
  private:
    std::shared_ptr<const SweeperInterface<F>> _ptr;
};

template<class SWP, class F> class SweeperBase
    : public virtual SweeperInterface<F>
{
    using typename SweeperInterface<F>::PR;
    virtual SweeperInterface<F>* _clone() const override final;
    virtual PR _precision() const override final;
    virtual Bool _discard(const MultiIndex& a, const F& x) const override final;
    virtual Void _sweep(Expansion<F>& p, F& e) const override;
    virtual Void _sweep(Expansion<F>& p) const override;
};


template<class SWP, class F>
SweeperInterface<F>* SweeperBase<SWP,F>::_clone() const
{
    return new SWP(static_cast<const SWP&>(*this));
}

template<class SWP, class F>
auto SweeperBase<SWP,F>::_precision() const -> PR
{
    return static_cast<const SWP*>(this)->SWP::precision();
}

template<class SWP, class F>
Bool SweeperBase<SWP,F>::_discard(const MultiIndex& a, const F& x) const
{
    return static_cast<const SWP*>(this)->SWP::discard(a,x);
}

template<class SWP, class F>
Void SweeperBase<SWP,F>::_sweep(Expansion<F>& p, F& e) const
{
    typename Expansion<F>::ConstIterator end=p.end();
    typename Expansion<F>::ConstIterator adv=p.begin();
    typename Expansion<F>::Iterator curr=p.begin();
    F::set_rounding_upward();
    F te=0.0;
    while(adv!=end) {
        if(this->_discard(adv->key(),adv->data())) {
            te+=abs(adv->data());
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    e+=te;
    p.resize(curr-p.begin());
    F::set_rounding_to_nearest();
}

template<class SWP, class F>
Void SweeperBase<SWP,F>::_sweep(Expansion<F>& p) const
{
    typename Expansion<F>::ConstIterator end=p.end();
    typename Expansion<F>::ConstIterator adv=p.begin();
    typename Expansion<F>::Iterator curr=p.begin();
    while(adv!=end) {
        if(this->_discard(adv->key(),adv->data())) {
        } else {
            *curr=*adv;
            ++curr;
        }
        ++adv;
    }
    p.resize(curr-p.begin());
}

//! \brief A sweeper class which discards terms whose absolute value is smaller than a threshold.
template<class F> class ThresholdSweeper : public SweeperBase<ThresholdSweeper<F>,F> {
    typedef PrecisionType<F> PR;
    PR _coefficient_precision;
    F _sweep_threshold;
  public:
    ThresholdSweeper(PR precision, F sweep_threshold)
        : _coefficient_precision(precision), _sweep_threshold(sweep_threshold) { ARIADNE_ASSERT(sweep_threshold>=0); }
    inline PR precision() const { return _coefficient_precision; }
    inline F sweep_threshold() const { return _sweep_threshold; }
    inline Bool discard(const MultiIndex& a, const F& x) const { return abs(x) < this->_sweep_threshold; }
  private:
    virtual Void _write(OutputStream& os) const { os << "ThresholdSweeper( sweep_threshold="<<this->_sweep_threshold<<" )"; };
};

//! \brief A sweeper class which does not discard any terms at all.
template<class F> class TrivialSweeper : public SweeperBase<TrivialSweeper<F>,F> {
    typedef PrecisionType<F> PR;
    PR _coefficient_precision;
  public:
    TrivialSweeper(PR precision) : _coefficient_precision(precision) { }
    inline PR precision() const { return _coefficient_precision; }
    inline Bool discard(const MultiIndex& a, const F& x) const { return false; }
  private:
    virtual Void _sweep(Expansion<F>& p, F& e) const final { }
    virtual Void _sweep(Expansion<F>& p) const final { }
    virtual Void _write(OutputStream& os) const { os << "TrivialSweeper"; }
};

//! \brief A sweeper class which only discards the zero term.
template<class F> class NullSweeper : public SweeperBase<NullSweeper<F>,F> {
    typedef PrecisionType<F> PR;
    PR _coefficient_precision;
  public:
    NullSweeper(PR precision) : _coefficient_precision(precision) { }
    inline PR precision() const { return _coefficient_precision; }
    inline Bool discard(const MultiIndex& a, const F& x) const { return x==0.0; }
  private:
    virtual Void _write(OutputStream& os) const { os << "NullSweeper"; }
};

//! \brief A sweeper class which discards non-affine terms.
template<class F> class AffineSweeper : public SweeperBase<AffineSweeper<F>,F> {
    typedef PrecisionType<F> PR;
    PR _coefficient_precision;
  public:
    AffineSweeper(PR precision) : _coefficient_precision(precision) { }
    inline PR precision() const { return _coefficient_precision; }
    inline Bool discard(const MultiIndex& a, const F& x) const { return a.degree()>1; }
  private:
    virtual Void _write(OutputStream& os) const { os << "AffineSweeper"; }
};

//! \brief A sweeper class which discards terms whose total degree is above some threshold.
template<class F> class GradedSweeper : public SweeperBase<GradedSweeper<F>,F> {
    typedef PrecisionType<F> PR;
    PR _coefficient_precision;
  public:
    GradedSweeper(PR precision, DegreeType degree) : _coefficient_precision(precision), _degree(degree) { }
    DegreeType degree() const { return this->_degree; }
    inline Bool discard(const MultiIndex& a, const Float64& x) const { return a.degree()>this->_degree; }
    inline Bool discard(const MultiIndex& a, const FloatMP& x) const { return a.degree()>this->_degree; }
  private:
    virtual Void _write(OutputStream& os) const { os << "GradedSweeper( degree="<<this->_degree<<" )"; }
  private:
    Nat _degree;
};

template<> inline Sweeper<Float64>::Sweeper() : _ptr(new ThresholdSweeper<Float64>(Precision64(),std::numeric_limits<float>::epsilon())) { }
template<> inline Sweeper<FloatMP>::Sweeper() : _ptr(new ThresholdSweeper<FloatMP>(PrecisionMP(64),std::numeric_limits<double>::epsilon())) { }


} // namespace Ariadne

#endif // ARIADNE_SWEEPER_HPP
