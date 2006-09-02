/***************************************************************************
 *            interval.h
 *
 *  Wed 4 May 2005
 *  Copyright 2005  Alberto Casagrande, Pieter Collins
 *  Email casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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
 
/*! \file interval.h
 *  \brief Intervals of real number types (currently implemented using Boost).
 */
 
#ifndef _ARIADNE_INTERVAL_H
#define _ARIADNE_INTERVAL_H

#include <iostream>

#include <boost/numeric/interval.hpp>
#include <boost/numeric/interval/io.hpp>

#include "../declarations.h"

#include "../numeric/numerical_traits.h"
#include "../numeric/arithmetic.h"
 
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
      return is;
    }
  }
}


namespace Ariadne {
  namespace Numeric {
    /*!
     * \brief A templated class representing an interval of real numbers.
     * \ingroup Numeric
     * 
     * An interval of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a val of real numbers with endpoints of type \a R.
     * All operations on an interval must be guarenteed to return an interval contining the exact result.
     * If \a T supports exact evaluation of a function, then the exact evaluation must be used.
     * If \a T is dense in the reals, e.g. dyadic or rational, then any approximate operations may be given a maximum error of computation.
     *
     * Currently implemented as a wrapper around the boost::numeric::interval class template from the Boost C++ library.
     */
    //using namespace boost::numeric::interval_lib;

//typedef boost::numeric::interval< double, policies< rounded_math<double>,checking_base<double> > > dInterval;
    template<typename R>
    class Interval : public
        boost::numeric::interval<R, 
            boost::numeric::interval_lib::policies< 
                boost::numeric::interval_lib::rounded_math<R>,
                boost::numeric::interval_lib::checking_base<R> 
            > 
        >
    {
      typedef boost::numeric::interval<R, 
            boost::numeric::interval_lib::policies< 
                boost::numeric::interval_lib::rounded_math<R>,
                boost::numeric::interval_lib::checking_base<R> 
            > 
        > _boost_interval;

     public:
      /*! \brief Default constructer constructs empty interval. */
      Interval()
        : _boost_interval(1,0) { }
      /*! \brief Construct from lower and upper bounds. */
      template<typename RL,typename RU> Interval(const RL& l, const RU& u)
        : _boost_interval(R(l),R(u)) { }
       /*! \brief Construct an interval with possibly different real type. */
      template<typename RX> Interval(const Interval<RX>& ivl)
        : _boost_interval(R(ivl.lower()),R(ivl.upper())) { }
     /*! \brief Construct a one-point interval. */
      template<typename RX> Interval(const RX& x)
        : _boost_interval(R(x),R(x)) { }
      /*! \brief Construct from a boost interval. */
      Interval(const _boost_interval& ivl)
        : _boost_interval(ivl) { }
      /*! \brief Construct from lower and upper bounds. */
      Interval(const R& l, const R& u)
        : _boost_interval(l,u) { }
      /*! \brief Assignment operator. */
      Interval<R>& operator=(const R& x) {
        *this=Interval<R>(x,x);
        return *this;
      }
      
      /*! \brief Equality operator. */
      bool operator==(const Interval<R>& ivl) const { 
        return this->_boost_interval::operator==(ivl); }
      /*! \brief Inequality operator. */
      bool operator!=(const Interval<R>& ivl) const{ 
        return this->_boost_interval::operator!=(ivl); }

      /*! \brief Comparison operator. */
      bool operator==(const R& x) const { 
        return this->_boost_interval::operator==(x); }
      /*! \brief Inequality operator. */
      bool operator!=(const R& x) const { 
        return this->_boost_interval::operator!=(x); }

      /*! \brief Tests if the interval is empty . */
      bool empty() const { return this->lower()>this->upper(); }
      /*! \brief Tests if the interval contains \a r. */
      bool contains(const R& r) const { return this->lower()<=r && r<=this->upper(); }
      /*! \brief Tests if the interior of the interval contains \a r. */
      bool interior_contains(const R& r) const { return this->lower()<r && r<this->upper(); }
      
      /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
      R centre() const { return (this->lower()+this->upper())/2; }
      /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
      R radius() const { return (this->upper()-this->lower())/2; }
      /*! \brief The length of the interval, given by \f$b-a\f$. */
      R length() const { return this->upper()-this->lower(); }
      
#ifdef DOXYGEN
      /*! \brief The interval of possible minima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> min(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible maxima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> max(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible absolute values of \a x in \a ivl. */
      friend Interval<R> abs(const Interval<R>& ivl);
      
      /*! \brief Interval negation. */
      friend Interval<R> operator-(const Interval<R>& ivl);
      /*! \brief Interval addition. */
      friend Interval<R> operator+(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Interval subtraction. */
      friend Interval<R> operator-(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Interval multiplication. */
      friend Interval<R> operator*(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Interval division. */
      friend Interval<R> operator/(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Integer power. */
      friend template<typename N> Interval<R> pow(const Interval<R>& x, const N& n);

      /*! \brief Tests equality. */
      friend bool equal<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests disjointness. */
      friend bool disjoint<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests intersection of interiors. */
      friend bool interiors_intersect<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 is a subset of \a ivl2. */
      friend bool subset<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief Tests if \a ivl1 is a subset of the interior of \a ivl2. */
      friend bool inner_subset<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      
      /*! \brief The intersection of \a ivl1 and \a ivl2. */
      friend Interval<R> intersection<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The closure of the intersection of the interiors of \a ivl1 and \a ivl2. */
      friend Interval<R> regular_intersection<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The smallest interval containing \a ivl1 and \a ivl2. */
      friend Interval<R> hull<>(const Interval<R>& ivl1, const Interval<R>& ivl2);
      
#endif
      
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
    };
    
    /* A reference to an interval. */
    template <typename R>
    class IntervalReference {
     public:
      IntervalReference(R& l, R& u) : _lower(&l), _upper(&u) { }
      IntervalReference(Interval<R>& ivl)
        : _lower(&const_cast<R&>(ivl.lower())), _upper(&const_cast<R&>(ivl.upper())) { }
      void operator=(const Interval<R>& ivl) { *_lower=ivl.lower(); *_upper=ivl.upper(); }
      operator Interval<R> () const { return Interval<R>(*_lower,*_upper); }
      R lower() const { return *_lower; }
      R upper() const { return *_upper; }
      R centre() const { return (*_lower+*_upper)/2; }
     private:
      R* _lower; R* _upper;
    };
    
    /* Numerical traits for interval template class. */
    template<typename R> class numerical_traits< Interval<R> > {
     public:
      typedef field_tag algebraic_category;
      typedef Interval<R> field_extension_type;
    };
    
    template<typename R>
    inline
    Interval<R> 
    min(const Interval<R>& ivl1, const Interval<R>& ivl2) {
      //std::cerr << "Interval::min<Interval<" << name<R>() << ">>" << std::endl;
      //return Interval<R>(min(ivl1.lower(),ivl2.lower()),min(ivl1.upper(),ivl2.upper()));
      return boost::numeric::min(ivl1,ivl2);
    }
  
    template<typename R>
    inline
    Interval<R> 
    max(const Interval<R>& ivl1, const Interval<R>& ivl2) {
      //std::cerr << "Interval::max<Interval<" << name<R>() << ">>" << std::endl;
      //return Interval<R>(max(ivl1.lower(),ivl2.lower()),max(ivl1.upper(),ivl2.upper()));
      return boost::numeric::max(ivl1,ivl2);
    }
    
    template<typename R>
    inline
    Interval<R> 
    abs(const Interval<R>& ivl) {
      using namespace ::Ariadne::Numeric;
      //std::cerr << "Interval::abs<Interval<" << name<R>() << ">>" << std::endl;
      //Can't use boost::numeric::abs because no comparison of GMP expressions.
      //return boost::numeric::abs(ivl);
      if(ivl.lower()>=0) { return ivl; } else if(ivl.upper() < 0) { return -ivl; } 
      else { return Interval<R>(0,max(R(-ivl.lower()),ivl.upper())); }
      //else { return Interval<R>(0,max<R>(-ivl.lower(),ivl.upper())); }
    }
  
    template<typename R>
    inline
    bool
    equal(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.lower()==ivl2.lower() || ivl1.upper()==ivl2.upper());
    }

    template<typename R>
    inline
    bool
    disjoint(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.upper()<ivl2.lower() || ivl1.lower()>ivl2.upper());
    }

    template<typename R>
    inline
    bool
    interiors_intersect(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.upper()>ivl2.lower() && ivl1.lower()<ivl2.upper());
    }

    template<typename R>
    inline
    bool
    inner_subset(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.lower()>ivl2.lower() && ivl1.upper()<ivl2.upper());
    }

    template<typename R>
    inline
    bool
    subset(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return (ivl1.lower()>=ivl2.lower() && ivl1.upper()<=ivl2.upper());
    }

    template<typename R>
    inline
    Interval<R>
    intersection(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      return Interval<R>(std::max(ivl1.lower(),ivl2.lower()),
                         std::min(ivl1.upper(),ivl2.upper()));
    }

    template<typename R>
    inline
    Interval<R>
    regular_intersection(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      if(ivl1.upper()==ivl2.lower() || ivl1.lower()==ivl2.upper()) {
        return Interval<R>(1,0);
      }
      else {
        return intersection(ivl1,ivl2);
      }
    }

    template<typename R>
    inline
    Interval<R>
    hull(const Interval<R>& ivl1, const Interval<R>& ivl2)
    {
      if(ivl1.empty()) {
        return ivl2;
      }
      if(ivl2.empty()) {
        return ivl1;
      }
      return Interval<R>(std::min(ivl1.lower(),ivl2.lower()),
                         std::max(ivl1.upper(),ivl2.upper()));
    }
    
    
    template<typename R,typename N>
    inline
    Interval<R> 
    pow(const Interval<R>& x, const N& n) 
    {
      Interval<R> result=R(1);
      for(N i=0; i!=n; ++i) {
        result*=x;
      }
      return result;
    }
  
    template<typename R>
    inline
    std::ostream&
    Interval<R>::write(std::ostream& os) const
    {
      if(this->empty()) {
        return os << "[1,0]";
      }
      else {
        return os << "[" << this->lower() << "," << this->upper() << "]";
      }
      //const _boost_interval& boost_ivl(*this);
      //return os << boost_ivl;
    }
    
    template<typename R>
    inline
    std::istream&
    Interval<R>::read(std::istream& is)
    {
      char c;
      R l;
      R u;
      is >> c >> l >> c >> u >> c;
      (*this)=Interval<R>(l,u);
      //_boost_interval& boost_ivl(*this);
      //is >> boost_ivl;
      return is;
    }
        
    template<typename R>
    inline
    std::ostream&
    operator<<(std::ostream& os, const Interval<R>& ivl)
    {
      return ivl.write(os);
    }
    
    template<typename R>
    inline
    std::istream&
    operator>>(std::istream& is, Interval<R>& ivl)
    {
      return ivl.read(is);
    }
    
  } // namespace Numeric
} // namespace Ariadne
  
#endif /* _ARIADNE_INTERVAL_H */
