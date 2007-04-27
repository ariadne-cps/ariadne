/***************************************************************************
 *            differential.h
 *
 *  Copyright 2007  Alberto Casagrande, Pieter Collins
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
 
/*! \file differential.h
 *  \brief Differentials for automatic differentiation.
 */
 
#ifndef ARIADNE_DIFFERENTIAL_H
#define ARIADNE_DIFFERENTIAL_H

#include <iostream>
#include <stdexcept>

#include "../base/tribool.h"
#include "../base/exceptions.h"

#include "../numeric/exceptions.h"
#include "../numeric/numerical_traits.h"
#include "../numeric/conversion.h"
#include "../numeric/arithmetic.h"
#include "../numeric/function.h"

namespace Ariadne {
  namespace Numeric {


    /*!\ingroup Numeric
     * \brief A templated class representing a the differential of a quantity.
     */
    template<class X, class V>
    class Differential
    {
     private:
      X _x; V _dx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs the differential of the constant zero with respect to no independent variable. */
      Differential() : _x(), _dx() { }
      /*! \brief Constuct a differential based on variable x and derivative vector dx. */
      template<class XX,class VV> Differential(const XX& x, const VV& dx) : _x(x), _dx(dx) { }
      /*! \brief Construct the differential of a constant with respect to \a n independent variables. */
      template<class XX> Differential(const XX& x, size_type n) : _x(x), _dx(n) { }
      /*! \brief Construct the differential of the \a i th variable of \a n independent variables. */
      template<class XX> Differential(const XX& x, size_type n, size_type i) : _x(x), _dx(n) { _dx[i]=1; }
      /*! \brief Assign the differential of a constant zero with respect to the current independent variables (useful as an alternative to a constructor). */
      template<class XX> Differential<X,V> operator=(const XX& c) { 
        this->_x=c; for(size_type j=0; j!=this->_dx.size(); ++j) { this->_dx[j]=0; } return *this; }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The differential of the variable with respect to the indepenent variables. */
      const V& derivative() const { return this->_dx; }
      /*! \brief The differential of the variable with respect to the \a j th indepenent variable. */
      const typename V::value_type& derivative(size_type j) const { return this->_dx[j]; }
      //@}
      
    };
    
    
    /*!\ingroup Numeric
     * \brief A templated class representing a the differential of a quantity with respect to a scalar argument.
     */
    template<class X>
    class Differential< X, X >
    {
     private:
      X _x; X _dx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor constructs a differential based on the constant zero. */
      Differential() : _x(), _dx() { }
      /*! \brief Constuct a differential based on variable x and derivative dx. */
      template<class XX0,class XX1> Differential(const XX0& x, const XX1& dx) : _x(x), _dx(dx) { }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The derivative of the variable with respect to the independent variable. */
      const X& derivative() const { return this->_dx; }
      //@}
    };
    


    template<class X1, class V1, class X2, class V2> inline
    bool operator==(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2) {
      return x1.value()==x2.value() && x1.derivative()==x2.derivative();
    }

    template<class X1, class V1, class X2, class V2> inline
    bool operator!=(const Differential<X1,V1>& x1, const Differential<X2,V2>& x2) {
      return !(x1==x2);
    }

    template<class X, class V> inline
    std::ostream& operator<<(std::ostream& os, const Differential<X,V>& x) {
      return os << "("<<x.value()<<","<<x.derivative()<<")";
    }


    template<class X> inline
    std::ostream& operator<<(std::ostream& os, const Differential< X, Differential<X,X> >& x) {
      return os << "("<<x.value()<<","<<x.first_derivative()<<","<<x.second_derivative()<<")";
    }



    
    template<class C, class X, class V> inline 
    Differential<X,V> operator*(const C& c, const Differential<X,V>& x) {
      return Differential<X,V>(c*x.value(),c*x.derivative());
    }
    
    template<class C, class X, class V> inline 
    Differential<X,V> operator*(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(c*x.value(),c*x.derivative());
    }
    
    template<class C, class X, class V> inline 
    Differential<X,V> operator/(const Differential<X,V>& x, const C& c) {
      return Differential<X,V>(x.value()/c,x.derivative()/c);
    }
    


    template<class X, class V> inline 
    Differential<X,V> min(const Differential<X,V>& x1,const Differential<X,V>& x2) {
      if(x1.value()==x2.value()) {
        ARIADNE_THROW(std::runtime_error,"min(Differntial x1, Differential x2)","x1.value()==x2.value()");
      }
      return x1.value()<x2.value() ? x1 : x2;
    }
    
    template<class X, class V> inline 
    Differential<X,V> max(const Differential<X,V>& x1,const Differential<X,V>& x2) {
      if(x1.value()==x2.value()) { 
        ARIADNE_THROW(std::runtime_error,"max(Differntial x1, Differential x2)","x1.value()==x2.value()"); 
      }
      return x1.value()>x2.value() ? x1 : x2;
    }
    
    template<class X, class V> inline 
    Differential<X,V> abs(const Differential<X,V>& x) {
      if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
      if(x.value()>0) { return x; } else { return -x; }
    }
    
    template<class X, class V> inline 
    Differential<X,V> pos(const Differential<X,V>& x) {
      return Differential<X,V>(x.value(),x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> neg(const Differential<X,V>& x) {
      return Differential<X,V>(-x.value(),-x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> add(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()+x2.value(),x1.derivative()+x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> sub(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()-x2.value(),x1.derivative()-x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> mul(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return Differential<X,V>(x1.value()*x2.value(),x1.derivative()*x2.value()+x1.value()*x2.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> div(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      X y=x1.value()/x2.value();
      V v=x1.derivative()-y*x2.derivative();
      return Differential<X,V>(y,v/x2.value());
    }
    
    template<class X, class V> inline 
    Differential<X,V> pow(const Differential<X,V>& x, int n) {
      //std::cerr << "pow(Differential x, int n)"<<std::endl;
      if(n==0) {
        return Differential<X,V>(pow(x.value(),0),static_cast<X>(0)*x.derivative());
      } else {
        X y=pow(x.value(),n-1);
        return Differential<X,V>((x.value()*y),static_cast<X>(n*y)*x.derivative());
      }
    }
    
    template<class X, class V> inline 
    Differential<X,V> sqrt(const Differential<X,V>& x) {
      X y=sqrt(x.value());
      return Differential<X,V>(y,x.derivative()/(2*y));
    }
    
    template<class X, class V> inline 
    Differential<X,V> exp(const Differential<X,V>& x) {
      X y=exp(x.value());
      return Differential<X,V>(y,y*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> log(const Differential<X,V>& x) {
      return Differential<X,V>(log(x.value()),x.derivative()/x.value());
    }
    
    template<class X, class V> inline 
    Differential<X,V> sin(const Differential<X,V>& x) {
      return Differential<X,V>(sin(x.value()),cos(x.value())*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> cos(const Differential<X,V>& x) {
      return Differential<X,V>(cos(x.value()),-sin(x.value())*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> tan(const Differential<X,V>& x) {
      X y=tan(x.value());
      return Differential<X,V>(y,(X(1)-y*y)*x.derivative());
    }
    
    template<class X, class V> inline 
    Differential<X,V> asin(const Differential<X,V>& x) {
      return Differential<X,V>(asin(x.value()),x.derivative()/sqrt(X(1)-x.value()*x.value()));
    }
    
    template<class X, class V> inline 
    Differential<X,V> acos(const Differential<X,V>& x) {
      return Differential<X,V>(acos(x.value()),-x.derivative()/sqrt(X(1)-x.value()*x.value()));
    }
    
    template<class X, class V> inline 
    Differential<X,V> atan(const Differential<X,V>& x) {
      return Differential<X,V>(atan(x.value()),x.derivative()/(X(1)+x.value()*x.value()));
    }
 


    template<class X, class V> inline 
    Differential<X,V> operator+(const Differential<X,V>& x) {
      return pos(x);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator-(const Differential<X,V>& x) {
      return neg(x);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator+(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return add(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator-(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return sub(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator*(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return mul(x1,x2);
    }
 
    template<class X, class V> inline 
    Differential<X,V> operator/(const Differential<X,V>& x1, const Differential<X,V>& x2) {
      return div(x1,x2);
    }







    template<class X,class V, class H> class SecondDifferential;

    /*!\ingroup Numeric
     * \brief A templated class representing a the second differential of a quantity with respect to a scalar argument.
     */
    template<class X>
    class SecondDifferential< X, X, X >
    {
      private:
      X _x; X _dx; X _ddx;
     public:
      //@{
      //! \name Constructors and assignment operators
      /*! \brief Default constructor. */
      SecondDifferential() : _x(), _dx(), _ddx() { }
      /*! \brief Constuct a second differential based on variable \a x, derivative \a dx and second derivative \a ddx. */
      template<class XX0,class XX1,class XX2> SecondDifferential(const XX0& x, const XX1& dx, const XX2& ddx)
        : _x(x), _dx(dx), _ddx(ddx) { }
      //@}

      //@{
      //! \name Data access
      /*! \brief The value of the variable. */
      const X& value() const { return this->_x; }
      /*! \brief The derivative of the variable with respect to the indepenent variable. */
      const X& derivative() const { return this->_dx; }
      /*! \brief The first derivative of the variable with respect to the indepenent variable. */
      const X& first_derivative() const { return this->_dx; }
      /*! \brief The second derivative of the variable with respect to the indepenent variable. */
      const X& second_derivative() const { return this->_ddx; }
      //@}
      
    };
    
      
      
    template<class X1, class V1, class H1, class X2, class V2, class H2> 
    bool operator==(const SecondDifferential<X1,V1,H1>& x1, const SecondDifferential<X2,V2,H2>& x2) {
      return x1.value()==x2.value() 
        && x1.derivative()==x2.derivative()
        && x1.second_derivative()==x2.second_derivative();
    }

    template<class X1, class V1, class H1, class X2, class V2, class H2> 
    bool operator!=(const SecondDifferential<X1,V1,H1>& x1, const SecondDifferential<X2,V2,H2>& x2) {
      return !(x1==x2);
    }

    template<class X, class V, class H> 
    std::ostream& operator<<(std::ostream& os, const SecondDifferential<X,V,H>& x) {
    return os << "("<<x.value()<<","<<x.derivative()<<","<<x.second_derivative()<<")";
    }




    template<class X>  
    SecondDifferential<X,X,X> abs(const SecondDifferential<X,X,X>& x) {
      if(x.value()==0) { ARIADNE_THROW(std::runtime_error,"abs(Differntial x)","x.value()==0"); }
      if(x.value()>0) { return pos(x); } else { return neg(x); }
    }
    
    template<class X>  
    SecondDifferential<X,X,X> pos(const SecondDifferential<X,X,X>& x) {
      return SecondDifferential<X,X,X>(x.value(),x.derivative(),x.second_derivative());
    }
    
    template<class X>  
    SecondDifferential<X,X,X> neg(const SecondDifferential<X,X,X>& x) {
      return SecondDifferential<X,X,X>(-x.value(),-x.derivative(),-x.second_derivative());
    }
    
    template<class X>  
    SecondDifferential<X,X,X> add(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()+x2.value();
      X dy=x1.derivative()+x2.derivative();
      X ddy=x1.second_derivative()+x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sub(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()-x2.value();
      X dy=x1.derivative()-x2.derivative();
      X ddy=x1.second_derivative()-x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> mul(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()*x2.value();
      X dy=x1.derivative()*x2.value()+x1.value()*x2.derivative();
      X ddy=x1.second_derivative()*x2.value()+2*x1.derivative()*x2.derivative()+x1.value()*x2.second_derivative();
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    // (ddx1*x2-2*x2*dx1*dx2-x1*ddx2)/x2^2+2*x1*dx2^2/x2^3
    template<class X>  
    SecondDifferential<X,X,X> div(const SecondDifferential<X,X,X>& x1, const SecondDifferential<X,X,X>& x2) {
      X y=x1.value()/x2.value();
      X s=x1.derivative()-y*x2.derivative();
      X dy=s/x2.value();
      X ds=x1.second_derivative()-dy*x2.derivative()-y*x2.second_derivative();
      X ddy=(ds-dy*x2.derivative())/x2.value();
      //std::cerr << "y="<<y<<" s="<<s<<" dy="<<dy<<" ds="<<ds<<" ddy="<<ddy<<std::endl;
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> pow(const SecondDifferential<X,X,X>& x, int n) {
      //std::cerr << "pow(SecondDifferential x, int n)"<<std::endl;
      if(n==0) {
        return SecondDifferential<X,X,X>(1,0,0);
      } else if(n==1) {
        return SecondDifferential<X,X,X>(x.value(),1,0);
      } else {
        X z2=(n*(n-1))*pow(x.value(),n-2);
        X z1=n*pow(x.value(),n-1);
        X z0=pow(x.value(),n);
        X w=pow(x.derivative(),2);
        return SecondDifferential<X,X,X>(z0,z1*x.derivative(),z1*x.second_derivative()+z2*w);
      }
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sqrt(const SecondDifferential<X,X,X>& x) {
      X y=sqrt(x.value());
      X r=static_cast<X>(0.5)/y;
      X dy=x.derivative()*r;
      X ddy=(x.second_derivative()-pow(x.derivative(),2)/(2*x.value()))*r;
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> exp(const SecondDifferential<X,X,X>& x) {
      X y=exp(x.value());
      X dy=y*x.derivative();
      X ddy=y*(pow(x.derivative(),2)+x.second_derivative());
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> log(const SecondDifferential<X,X,X>& x) {
      X y=log(x.value());
      X dy=x.derivative()/x.value();
      X ddy=x.second_derivative()/x.value()-pow(dy,2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> sin(const SecondDifferential<X,X,X>& x) {
      X s=sin(x.value());
      X c=cos(x.value());
      const X& y=s;
      X dy=c*x.derivative();
      X ddy=c*x.second_derivative()-s*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> cos(const SecondDifferential<X,X,X>& x) {
      X c=cos(x.value());
      X s=sin(x.value());
      const X& y=c;
      X dy=-s*x.derivative();
      X ddy=-s*x.second_derivative()-c*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> tan(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X y=tan(x.value());
      X z=1+pow(y,2);
      X dy=z*x.derivative();
      X ddy=dy*x.second_derivative()+2*y*z*pow(x.derivative(),2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }
    
    template<class X>  
    SecondDifferential<X,X,X> asin(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      // TODO: Check these formulae
      X r=1/sqrt(1-pow(x.value(),2));
      X y=asin(x.value());
      X dy=r*x.derivative();
      X ddy=r*(x.second_derivative()+pow(dy,2));
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }


    template<class X>  
    SecondDifferential<X,X,X> acos(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X r=1/sqrt(1-pow(x.value(),2));
      X y=acos(x.value());
      X dy=r*x.derivative();
      X ddy=r*(x.second_derivative()+pow(dy,2));
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }

    template<class X>  
    SecondDifferential<X,X,X> atan(const SecondDifferential<X,X,X>& x) {
      // TODO: Check these formulae
      X r=1/(1+pow(x.value(),2));
      X y=atan(x.value());
      X dy=x.derivative()*r;
      X ddy=x.second_derivative()*r-2*x.value()*pow(dy,2);
      return SecondDifferential<X,X,X>(y,dy,ddy);
    }


    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator+(const SecondDifferential<X,V,H>& x) {
      return pos(x);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator-(const SecondDifferential<X,V,H>& x) {
      return neg(x);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator+(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return add(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator-(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return sub(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator*(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return mul(x1,x2);
    }
 
    template<class X, class V, class H> inline 
    SecondDifferential<X,V,H> operator/(const SecondDifferential<X,V,H>& x1, const SecondDifferential<X,V,H>& x2) {
      return div(x1,x2);
    }



  }
}

#endif /* ARIADNE_DIFFERENTIAL_H */

 
