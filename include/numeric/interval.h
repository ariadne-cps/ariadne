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
#include <stdexcept>

#include "../declarations.h"

#include "../base/tribool.h"

#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Numeric {
    //using namespace boost::numeric::interval_lib;

  class DivideByZeroException : public std::exception { };

    /*!\ingroup Numeric
     * \brief A templated class representing an interval of real numbers.
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
    template<class R>
    class Interval
    {
     private:
      R _lower; R _upper;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructer constructs empty interval. */
      Interval() : _lower(conv_down<R>(0)), _upper(conv_up<R>(0)) { }
      /*! \brief Construct from lower and upper bounds. */
      template<class RL,class RU> Interval(const RL& l, const RU& u)
        : _lower(conv_down<R>(l)), _upper(conv_up<R>(u)) { }
       /*! \brief Construct an interval with possibly different real type. */
      template<class RX> Interval(const Interval<RX>& ivl)
        : _lower(conv_down<R>(ivl.lower())), _upper(conv_up<R>(ivl.upper())) { }
      /*! \brief Construct a one-point interval. */
      template<class RX> Interval(const RX& x)
        : _lower(conv_down<R>(x)), _upper(conv_up<R>(x)) { }
      /*! \brief Construct from lower and upper bounds. */
      Interval(const R& l, const R& u) : _lower(l), _upper(u) { }
      
      /*! \brief Assignment operator. */
      template<class RX> 
      Interval<R>& operator=(const RX& x) {
        this->_lower=conv_down<R>(x); this->_upper=conv_up<R>(x); return *this;
      }
      /*! \brief Assignment operator. */
      Interval<R>& operator=(const R& x) {
        this->_lower=x; this->_upper=x; return *this;
      }
      /*! \brief Assignment operator. */
      template<class RX> 
      Interval<R>& operator=(const Interval<RX>& ivl) {
        this->_lower=conv_down<R>(ivl.lower()); this->_upper=conv_up<R>(ivl.upper()); return *this;
      }
      /*! \brief Copy assignment operator. */
      Interval<R>& operator=(const Interval<R>& ivl) {
        this->_lower=ivl._lower; this->_upper=ivl._upper; return *this;
      }
      //@}
      
      //@{
      //! \name Data access
      /*! \brief The lower bound. */
      const R& lower() const { return this->_lower; }
      /*! \brief The upper bound. */
      const R& upper() const { return this->_upper; }
      //@}
      
      //@{
      //! \name Geometric operations
      /*! \brief The midpoint of the interval, given by \f$(a+b)/2\f$. */
      R centre() const { 
        return div_approx(add_approx(this->lower(),this->upper()),R(2)); }
      /*! \brief The radius of the interval, given by \f$(b-a)/2\f$. */
      R radius() const { return div_up(sub_up(this->upper(),this->lower()),R(2)); }
      /*! \brief The length of the interval, given by \f$b-a\f$. */
      R length() const { return sub_up(this->upper(),this->lower()); }
      
      /*! \brief Tests if the interval is empty. */
      bool empty() const { return this->lower()>this->upper(); }
      /*! \brief Tests if the interval consists of a single point. */
      bool singleton() const { return this->lower()==this->upper(); }
      /*! \brief Tests if the interval contains \a r. */
      bool contains(const R& r) const { return this->lower()<=r && r<=this->upper(); }
      //@}
      
#ifdef DOXYGEN
      //@{
      //! \name Arithmetic operations
      /*! \brief The interval of possible minima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> min(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible maxima of \a x1 in \a ivl1 and \a x2 in \a ivl2. */
      friend Interval<R> max(const Interval<R>& ivl1, const Interval<R>& ivl2);
      /*! \brief The interval of possible absolute values of \a x in \a ivl. */
      friend Interval<R> abs(const Interval<R>& ivl);
      
      /*! \brief In-place addition of an interval. */
      friend Interval<R>& operator+=<>(Interval<R>&, const Interval<R>&);
      /*! \brief In-place addition of a number. */
      friend Interval<R>& operator+=<>(Interval<R>&, const R&);
      /*! \brief In-place subtraction of an interval. */
      friend Interval<R>& operator-=<>(Interval<R>&, const Interval<R>&);
      /*! \brief In-place subtraction of a number. */
      friend Interval<R>& operator-=<>(Interval<R>&, const R&);

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
      friend template<class N> Interval<R> pow(const Interval<R>& x, const N& n);
      //@}
      
      //@{
      //! \name Geometric operations
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
      //@}
#endif
      
      //@{
      //! \name Input/output operations
      /*! \brief Write to an output stream . */
      std::ostream& write(std::ostream& os) const;
      /*! \brief Read from an input stream . */
      std::istream& read(std::istream& is);
      //@}
    };
    
    /*!\ingroup Numeric
     * \brief A reference to an interval. 
     */
    template<class R>
    class IntervalReference {
     public:
      IntervalReference(R& l, R& u) : _lower(&l), _upper(&u) { }
      IntervalReference(Interval<R>& ivl)
        : _lower(&const_cast<R&>(ivl.lower())), _upper(&const_cast<R&>(ivl.upper())) { }
      void operator=(const Interval<R>& ivl) { *_lower=ivl.lower(); *_upper=ivl.upper(); }
      operator Interval<R> () const { return Interval<R>(*_lower,*_upper); }
      R lower() const { return *_lower; }
      R upper() const { return *_upper; }
      R centre() const { return Interval<R>(*this).centre(); }
     private:
      R* _lower; R* _upper;
    };
    
    template<class R> inline
    R approximate_value(const Interval<R>& ivl) 
    {
      return ivl.centre();
      //return Numeric::med_appox(ivl.lower(),ivl.upper());
      //return div_approx(add_appox(ivl.lower(),ivl.upper()),2);
    }
    
    template<class R> inline
    bool contains_value(const Interval<R>& ivl, const R& x) 
    {
      return ivl.contains(x);
    }
    
    template<class R> inline
    R error_bound(const Interval<R>& ivl) 
    {
      R av=approximate_value(ivl);
      return max_up(sub_up(av,ivl.lower()),sub_up(ivl.upper(),av));
    }
    



    /*! \brief Equality operator. */
    template<class R1, class R2> inline
    tribool operator==(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()==ivl2.upper() && ivl1.upper()==ivl2.lower()) { return true; }
      if(ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower()) { return false; }
      return indeterminate;
    }
      
    /*! \brief Inequality operator. */
    template<class R1, class R2> inline
    tribool operator!=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower()) { return true; }
      if(ivl1.lower()==ivl2.upper() && ivl1.upper()==ivl2.lower()) { return false; }
      return indeterminate;
    }
        
    /*! \brief Less than operator. */
    template<class R1, class R2> inline
    tribool operator<(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.upper()<ivl2.lower()) { return true; } 
      if(ivl1.lower()>=ivl2.upper()) { return false; }
      return indeterminate;
    }
        
    /*! \brief Greater than operator. */
    template<class R1, class R2> inline
    tribool operator>(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>ivl2.upper()) { return true; }
      if(ivl1.upper()<=ivl2.lower()) { return false; }
      return indeterminate;
    }
      
    /*! \brief Less than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator<=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.upper()<=ivl2.lower()) { return true; } 
      if(ivl1.lower()>ivl2.upper()) { return false; }
      return indeterminate;
    }
        
    /*! \brief Greater than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator>=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>=ivl2.upper()) { return true; }
      if(ivl1.upper()<ivl2.lower()) { return false; }
      return indeterminate;
    }
      


    /*! \brief Equality operator. */
    template<class R1, class R2> inline
    tribool operator==(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()==x() && ivl.upper()==x) { return true; }
      if(ivl.lower()>x || ivl.upper()<x) { return false; }
      return indeterminate;
    }
      
    /*! \brief inquality operator. */
    template<class R1, class R2> inline
    tribool operator!=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>x || ivl.upper()<x) { return true; }
      if(ivl.lower()==x() && ivl.upper()==x) { return false; }
      return indeterminate;
    }
      
    /*! \brief Less than operator. */
    template<class R1, class R2> inline
    tribool operator<(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.upper()<x) { return true; }
      if(ivl.lower()>=x) { return false; }
      return indeterminate;
    }
      
    /*! \brief Greater than operator. */
    template<class R1, class R2> inline
    tribool operator>(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>x) { return true; }
      if(ivl.upper()<=x) { return false; }
      return indeterminate;
    }
      
    /*! \brief Less than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator<=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()<=x) { return true; }
      if(ivl.upper()>x) { return false; }
      return indeterminate;
    }

    /*! \brief Greater than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator>=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>=x) { return true; }
      if(ivl.upper()<x) { return false; }
      return indeterminate;
    }


      
    /*! \brief Equality operator. */
    template<class R1, class R2> inline
    tribool operator==(const R1& x, const Interval<R2>& ivl) { 
      return ivl==x; 
    }

      
    /*! \brief inquality operator. */
    template<class R1, class R2> inline
    tribool operator!=(const R1& x, const Interval<R2>& ivl) {
      return ivl!=x;
    }
      
    /*! \brief Less than operator. */
    template<class R1, class R2> inline
    tribool operator<(const R1& x, const Interval<R2>& ivl) {
      return ivl>x;
    }
      
    /*! \brief Greater than operator. */
    template<class R1, class R2> inline
    tribool operator>(const R1& x, const Interval<R2>& ivl) {
      return ivl<x;
    }

    /*! \brief Less than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator<=(const R1& x, const Interval<R2>& ivl) {
      return ivl>=x;
    }
      
    /*! \brief Greater than or equal to operator. */
    template<class R1, class R2> inline
    tribool operator>=(const R1& x, const Interval<R2>& ivl) {
      return ivl<=x;
    }










    template<class R> inline
    Interval<R> operator+(const Interval<R>& x) {
      return x;
    }

    template<class R> inline
    Interval<R> operator-(const Interval<R>& x) {
      return Interval<R>(neg_down(x.upper()),neg_up(x.lower()));
    }


    
    template<class R> inline
    Interval<R> operator+(const Interval<R>& x1, const Interval<R>& x2) {
      //std::cerr << "operator+(" << name< Interval<R> >() << "," << name< Interval<R> >() << ")\n";
      return Interval<R>(add_down(x1.lower(),x2.lower()),
                         add_up(x1.upper(),x2.upper()));
    }

    template<class R> inline
    Interval<R> operator+(const Interval<R>& x1, const R& x2) {
      //std::cerr << "operator+(" << name< Interval<R> >() << "," << name<R>() << ")\n";
      return Interval<R>(add_down(x1.lower(),x2),add_up(x1.upper(),x2));
    }

    template<class R> inline
    Interval<R> operator+(const R& x1, const Interval<R>& x2) {
      //std::cerr << "operator+(" << name<R>() << "," << name< Interval<R> >() << ")\n";
      return Interval<R>(add_down(x1,x2.lower()),add_up(x1,x2.upper()));
    }


    template<class R> inline
    Interval<R> operator+=(Interval<R>& x1, const Interval<R>& x2) {
      //std::cerr << "operator+=(" << name< Interval<R> >() << "," << name< Interval<R> >() << ")\n";
      return x1=x1+x2;
    }

    template<class R> inline
    Interval<R> operator+=(Interval<R>& x1, const R& x2) {
      //std::cerr << "operator+=(" << name< Interval<R> >() << "," << name<R>() << ")\n";
      return x1=x1+x2;
    }



    template<class R> inline
    Interval<R> operator-(const Interval<R>& x1, const Interval<R>& x2) {
      return Interval<R>(sub_down(x1.lower(),x2.upper()),
                         sub_up(x1.upper(),x2.lower()));
    }

    template<class R> inline
    Interval<R> operator-(const Interval<R>& x1, const R& x2) {
      return Interval<R>(sub_down(x1.lower(),x2),sub_up(x1.upper(),x2));
    }

    template<class R> inline
    Interval<R> operator-(const R& x1, const Interval<R>& x2) {
      return Interval<R>(sub_down(x1,x2.upper()),sub_up(x1,x2.lower()));
    }

    template<class R> inline
    Interval<R> operator-=(Interval<R>& x1, const Interval<R>& x2) {
      return x1=x1-x2;
    }

    template<class R> inline
    Interval<R> operator-=(Interval<R>& x1, const R& x2) {
      return x1=x1-x2;
    }



    template<class R> inline
    Interval<R> operator*(const Interval<R>& x1, const Interval<R>& x2) {
      typedef Interval<R> I;
      const R& xl = x1.lower();
      const R& xu = x1.upper();
      const R& yl = x2.lower();
      const R& yu = x2.upper();
      R z=0;
      if (xl>=z) {
        if (yl>=z) {
          return I(mul_down(xl,yl),mul_up(xu,yu));
        } else if(yu<=z) {
          return I(mul_down(xu,yl),mul_up(xl,yu));
        } else {
          return I(mul_down(xu,yl),mul_up(xu,yu));
        }
      } else if (xu<=z) {
        if (yl>=z) {
          return I(mul_down(xl,yu),mul_up(xu,yl));
        } else if(yu<=z) {
          return I(mul_down(xu,yu),mul_up(xl,yl));
        } else {
          return I(mul_down(xl,yu),mul_up(xl,yl));
        }
      } else {
        if (yl>=z) {
          return I(mul_down(xl,yu),mul_up(xu,yu));
        } else if(yu<=z) {
          return I(mul_down(xu,yl),mul_up(xl,yl));
        } else {
          return I(min_down(mul_down(xl,yu),mul_down(xu,yl)),
                   max_up(mul_up(xl,yl),mul_up(xu,xu)));
        }
      }
    }    

    template<class R> inline
    Interval<R> operator*(const Interval<R>& x1, const R& x2) {
      return x1*Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator*(const R& x1, const Interval<R>& x2) {
      return Interval<R>(x1)*x2;
    }

    template<class R> inline
    Interval<R> operator*=(Interval<R>& x1, const Interval<R>& x2) {
      return x1=x1*x2;
    }

    template<class R> inline
    Interval<R> operator*=(Interval<R>& x1, const R& x2) {
      return x1=x1*x2;
    }


    template<class R> inline
    Interval<R> operator/(const Interval<R>& x1, const Interval<R>& x2) {
      typedef Interval<R> I;
      const R& xl = x1.lower();
      const R& xu = x1.upper();
      const R& yl = x2.lower();
      const R& yu = x2.upper();
      R z=0;
      if (yl>z) {
        if (xl>=z) {
          return I(div_down(xl,yu),div_up(xu,yl));
        } else if(xu<=z) {
          return I(div_down(xl,yl),div_up(xu,yu));
        } else {
          return I(div_down(xl,yl),div_up(xu,yl));
        }
      } else if (yu<z) {
        if (xl>=z) {
          return I(div_down(xu,yu),div_up(xl,yl));
        } else if(xu<=z) {
          return I(div_down(xu,yl),div_up(xl,yu));
        } else {
          return I(div_down(xu,yu),div_up(xl,yu));
        }
      } else {
        throw DivideByZeroException();
      }
    }

    template<class R> inline
    Interval<R> operator/(const Interval<R>& x1, const R& x2) {
      return x1/Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator/(const R& x1, const Interval<R>& x2) {
      return Interval<R>(x1)/x2;
    }

    template<class R> inline
    Interval<R> operator/=(Interval<R>& x1, const Interval<R>& x2) {
      return x1=x1/x2;
    }

    template<class R> inline
    Interval<R> operator/=(Interval<R>& x1, const R& x2) {
      return x1=x1/x2;
    }




    template<class R> inline
    Interval<R> min(const Interval<R>& x1, const Interval<R>& x2) {
      //std::cerr << "Interval::min<Interval<" << name<R>() << ">>" << std::endl;
      return Interval<R>(min_down(x1.lower(),x2.lower()),min_up(x1.upper(),x2.upper()));
    }
    
    template<class R> inline
    Interval<R> max(const Interval<R>& x1, const Interval<R>& x2) {
      //std::cerr << "Interval::max<Interval<" << name<R>() << ">>" << std::endl;
      return Interval<R>(max_down(x1.lower(),x2.lower()),max_up(x1.upper(),x2.upper()));
    }

    template<class R> inline
    Interval<R> abs(const Interval<R>& x) {
      using namespace ::Ariadne::Numeric;
      if(x.lower()>=0) { return x; } 
      if(x.upper() < 0) { return -x; } 
      return Interval<R>(R(0),max_up(neg_up(x.lower()),x.upper()));
    }
  
    template<class R,class N> inline
    Interval<R> pow(const Interval<R>& x, const N& n) {
      Interval<R> result=R(1);
      for(N i=0; i!=n; ++i) {
        result*=x;
      }
      return result;
    }
  


    template<class R> inline
    R centre(const Interval<R>& x) {
      return x.centre();
    }
    
    template<class R> inline
    R radius(const Interval<R>& x) { 
      return x.radius();
    }
      
    template<class R> inline
    R length(const Interval<R>& x) { 
      return x.length();
    }


    template<class R> inline
    bool equal(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.lower()==x2.lower() && x1.upper()==x2.upper());
    }

    template<class R> inline
    bool disjoint(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.upper()<x2.lower() || x1.lower()>x2.upper());
    }

    template<class R> inline
    bool subset(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.lower()>=x2.lower() && x1.upper()<=x2.upper());
    }

    template<class R> inline
    Interval<R> intersection(const Interval<R>& x1, const Interval<R>& x2) {
      return Interval<R>(max_down(x1.lower(),x2.lower()),
                         min_up(x1.upper(),x2.upper()));
    }

    template<class R> inline
    Interval<R> hull(const Interval<R>& x1, const Interval<R>& x2) {
      if(x1.empty()) { return x2; }
      if(x2.empty()) { return x1; }
      return Interval<R>(min_down(x1.lower(),x2.lower()),
                         max_up(x1.upper(),x2.upper()));
    }
    
   

    
    


    template<class R> inline
    Interval<R> sqrt(const Interval<R>& x) {
      return Interval<R>(sqrt_down(x.lower()),sqrt_up(x.upper()));
    }

    template<class R> inline
    Interval<R> hypot(const Interval<R>& x1, const Interval<R>& x2) {
      return Interval<R>(hypot_down(x1.lower(),x2.lower()),
                         hypot_up(x1.upper(),x2.upper()));
    }

    template<class R> inline
    Interval<R> exp(const Interval<R>& x) {
      return Interval<R>(exp_down(x.lower()),exp_up(x.upper()));
    }

    template<class R> inline
    Interval<R> log(const Interval<R>& x) {
      return Interval<R>(log_down(x.lower()),log_up(x.upper()));
    }

    
    template<typename R> inline 
    Interval<R> sin(const Interval<R>& ivl) {
      R pl=pi_down<R>();
      R pu=pi_up<R>();
      int n=quot(ivl.lower(),pu);
      R l=sub_down(ivl.lower(),mul_up(pu,n));
      R u=sub_up(ivl.upper(),mul_down(pl,n));
      if(sub_down(u,l)>=mul_up(2,pu)) {
        return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
      } else {
        throw std::runtime_error("sin(const Interval<R>& ivl) not implemented completely");
      }
    }
    
    
    
    template<class R> inline
    std::ostream& Interval<R>::write(std::ostream& os) const {
      if(this->empty()) {
        return os << "[1,0]";
      }
      else {
        return os << "[" << this->lower() << "," << this->upper() << "]";
      }
    }
    
    
    template<class R> inline
    std::istream& Interval<R>::read(std::istream& is) {
      char c;
      R l;
      R u;
      is >> c >> l >> c >> u >> c;
      (*this)=Interval<R>(l,u);
      return is;
    }
        
    
    
    template<class R> inline
    std::ostream& operator<<(std::ostream& os, const Interval<R>& x) {
      return x.write(os);
    }
    
    
    template<class R> inline
    std::istream& operator>>(std::istream& is, Interval<R>& x) {
      return x.read(is);
    }
    
  } // namespace Numeric
} // namespace Ariadne
  

namespace TBLAS {
  template<class real>
  int
  iamax (const int N, const real *X, const int incX);

  template<class real>
  int
  iamax (const int N, const Ariadne::Numeric::Interval<real> *X, const int incX)
  {
#ifdef DEBUG
    std::cerr << "TBLAS::iamax(const int N, const interval<real> *X, const int incX)\n";
#endif

    real mx = 0;
    int ix = 0;
    int i;
    int result = 0;
    
    if (incX <= 0) {
      return 0;
    }

    for (i = 0; i < N; i++) {
      real av=Ariadne::Numeric::min_down(
                Ariadne::Numeric::abs_down(X[ix].lower()),
                Ariadne::Numeric::abs_down(X[ix].upper()));
      if (av > mx) {
        mx = av;
        result = i;
      }
      ix += incX;
    }
    
    return result;
  }
  
}


#endif /* _ARIADNE_INTERVAL_H */
