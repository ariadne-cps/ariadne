/***************************************************************************
 *            interval.inline.h
 *
 *  Copyright 2005-7  Alberto Casagrande, Pieter Collins
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

#include "numeric/expression.h"
#include "numeric/rounding.h"
#include "numeric/transcendental.h"
 
namespace Ariadne {
  namespace Numeric {
  
    using std::min;
    using std::max;


    template<class R> inline 
    Interval<R>::Interval()
      : _lower(0), _upper(0) { }

    template<class R> inline 
    Interval<R>::Interval(const R& x)
      : _lower(x), _upper(x) { }

    template<class R> inline 
    Interval<R>::Interval(const R& l, const R& u)
      : _lower(l), _upper(u) { }

    template<class R> inline 
    Interval<R>::Interval(const Interval<R>& ivl)
      : _lower(ivl._lower), _upper(ivl._upper) { }
  
    template<class R> inline 
    Interval<R>& Interval<R>::operator=(const R& x) {
      this->_lower=x; this->_upper=x; return *this;
    }

    template<class R> inline 
    Interval<R>& Interval<R>::operator=(const Interval<R>& ivl) {
      this->_lower=ivl._lower; this->_upper=ivl._upper; return *this;
    }

    template<class R> template<class RL,class RU> inline 
    Interval<R>::Interval(const RL& l, const RU& u)
      : _lower(l,round_down), _upper(u,round_up) { }

    template<class R> template<class RX> inline 
    Interval<R>::Interval(const Interval<RX>& ivl)
      : _lower(ivl.lower(),round_down), _upper(ivl.upper(),round_up) { }

    template<class R> template<class RX> inline 
    Interval<R>::Interval(const RX& x)
      : _lower(x,round_down), _upper(x,round_up) { }
      
    template<class R> template<class RX> inline 
    Interval<R>& Interval<R>::operator=(const RX& x) {
      set_(this->_lower,x,round_down); set_(this->_upper,x,round_up); return *this;
    }
  

    template<class R> template<class RX> inline 
    Interval<R>& Interval<R>::operator=(const Interval<RX>& ivl) {
      set_(this->_lower,ivl.lower(),round_down); set_(this->_upper,ivl.upper(),round_up); return *this;
    }


      
    template<class R> inline 
    const R& Interval<R>::lower() const { 
      return this->_lower; 
    }

    template<class R> inline 
    const R& Interval<R>::upper() const { 
      return this->_upper; 
    }
      

    template<class R> inline 
    R Interval<R>::midpoint() const { 
      return div<RoundApprox>(add<RoundApprox>(this->lower(),this->upper()),2);
    }

    template<class R> inline 
    R Interval<R>::radius() const { 
      return div<RoundUp>(sub<RoundUp>(this->upper(),this->lower()),2); 
    }

    template<class R> inline
    R Interval<R>::width() const { 
      return sub<RoundUp>(this->upper(),this->lower()); 
    }
      

    template<class R> inline 
    bool Interval<R>::empty() const { 
      return this->lower()>this->upper(); 
    }

    template<class R> inline 
    bool Interval<R>::singleton() const { 
      return this->lower()==this->upper();
    }

    template<class R> template<class RX> inline 
    bool Interval<R>::encloses(const RX& x) const { 
      return this->lower()<=x && x<=this->upper();
    }
      
    template<class R> template<class RX> inline 
    bool Interval<R>::refines(const Interval<RX>& ivl) const { 
      return ivl.lower()<=this->lower() && this->upper()<=ivl.upper(); 
    }
      



    template<class R> inline 
    IntervalReference<R>::IntervalReference(R& l, R& u)
      : _lower(&l), _upper(&u) { }

    template<class R> inline 
    IntervalReference<R>::IntervalReference(Interval<R>& ivl)
      : _lower(&const_cast<R&>(ivl.lower())), _upper(&const_cast<R&>(ivl.upper())) { }

    template<class R> inline
    void IntervalReference<R>::operator=(const Interval<R>& ivl) {
      *_lower=ivl.lower(); *_upper=ivl.upper(); }

    template<class R> inline
    IntervalReference<R>::operator Interval<R> () const {
      return Interval<R>(*_lower,*_upper); }
      
    template<class R> inline 
    R IntervalReference<R>::lower() const { 
      return *_lower; }

    template<class R> inline 
    R IntervalReference<R>::upper() const { 
      return *_upper; }
    
    template<class R> inline 
    R IntervalReference<R>::centre() const { 
      return Interval<R>(*this).centre(); }
    
    



    template<class R> inline
    bool empty(const Interval<R>& x) {
      return x.empty();
    }
    
    template<class R> inline
    R lower(const Interval<R>& x) {
      return x.lower();
    }
    
    template<class R> inline
    R upper(const Interval<R>& x) {
      return x.upper();
    }
    
    template<class R> inline
    R midpoint(const Interval<R>& x) {
      return x.midpoint();
    }
    
    template<class R> inline
    R radius(const Interval<R>& x) { 
      return x.radius();
    }
      
    template<class R> inline
    R width(const Interval<R>& x) { 
      return x.width();
    }


    template<class R> inline
    R mignitude(const Interval<R>& x) { 
      if(x.lower()>0) { return x.lower(); }
      else if(x.upper()<0) { return -x.upper(); }
      else { return static_cast<R>(0); }
    }


    template<class X, class Y> inline
    bool encloses(const Interval<X>& x, const Interval<Y>& y) { 
      return x.lower()<=y.lower() && x.upper()>=y.upper();
    }

    template<class X, class Y> inline
    bool encloses(const Interval<X>& ivl, const Y& y) { 
      return ivl.encloses(y);
    }

    template<class X, class Y> inline
    bool refines(const Interval<X>& x, const Interval<Y>& y) { 
      return x.refines(y);
    }


    template<class X, class Y> inline
    tribool operator==(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.lower()==y.upper() && x.upper()==y.lower()) { return true; }
      if(x.lower()>y.upper() || x.upper()<y.lower()) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator!=(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.lower()>y.upper() || x.upper()<y.lower()) { return true; }
      if(x.lower()==y.upper() && x.upper()==y.lower()) { return false; }
      return indeterminate;
    }
        
    
    template<class X, class Y> inline
    tribool operator<(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.upper()<y.lower()) { return true; } 
      if(x.lower()>=y.upper()) { return false; }
      return indeterminate;
    }
        
    
    template<class X, class Y> inline
    tribool operator>(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.lower()>y.upper()) { return true; }
      if(x.upper()<=y.lower()) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator<=(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.upper()<=y.lower()) { return true; } 
      if(x.lower()>y.upper()) { return false; }
      return indeterminate;
    }
        
    
    template<class X, class Y> inline
    tribool operator>=(const Interval<X>& x, const Interval<Y>& y) { 
      if(x.lower()>=y.upper()) { return true; }
      if(x.upper()<y.lower()) { return false; }
      return indeterminate;
    }
      


    
    template<class X, class Y> inline
    tribool operator==(const Interval<X>& ivl, const Y& x) { 
      if(ivl.lower()==x && ivl.upper()==x) { return true; }
      if(ivl.lower()>x || ivl.upper()<x) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator!=(const Interval<X>& ivl, const Y& x) { 
      if(ivl.lower()>x || ivl.upper()<x) { return true; }
      if(ivl.lower()==x && ivl.upper()==x) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator<(const Interval<X>& ivl, const Y& x) { 
      if(ivl.upper()<x) { return true; }
      if(ivl.lower()>=x) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator>(const Interval<X>& ivl, const Y& x) { 
      if(ivl.lower()>x) { return true; }
      if(ivl.upper()<=x) { return false; }
      return indeterminate;
    }
      
    
    template<class X, class Y> inline
    tribool operator<=(const Interval<X>& ivl, const Y& x) { 
      if(ivl.upper()<=x) { return true; }
      if(ivl.lower()>x) { return false; }
      return indeterminate;
    }

    
    template<class X, class Y> inline
    tribool operator>=(const Interval<X>& ivl, const Y& x) { 
      if(ivl.lower()>=x) { return true; }
      if(ivl.upper()<x) { return false; }
      return indeterminate;
    }

    
    template<class X, class Y> inline
    tribool operator==(const X& x, const Interval<Y>& ivl) { 
      return ivl==x; 
    }

      
    
    template<class X, class Y> inline
    tribool operator!=(const X& x, const Interval<Y>& ivl) {
      return ivl!=x;
    }
      
    
    template<class X, class Y> inline
    tribool operator<(const X& x, const Interval<Y>& ivl) {
      return ivl>x;
    }
      
    
    template<class X, class Y> inline
    tribool operator>(const X& x, const Interval<Y>& ivl) {
      return ivl<x;
    }

    
    template<class X, class Y> inline
    tribool operator<=(const X& x, const Interval<Y>& ivl) {
      return ivl>=x;
    }
      
    
    template<class X, class Y> inline
    tribool operator>=(const X& x, const Interval<Y>& ivl) {
      return ivl<=x;
    }



    template<class R, class X> inline
    void pos_(Interval<R>& r, const Interval<X>& x) {
      pos_(r._lower,x._lower); 
      pos_(r._upper,x._upper);
    }

    template<class R, class X> inline
    void pos_(Interval<R>& r, const X& x) {
      pos_(r._lower,x); 
      pos_(r._upper,x);
    }


    template<class R, class X> inline
    void neg_(Interval<R>& r, const Interval<X>& x) {
      neg_(r._lower,x._upper);
      neg_(r._upper,x._lower);
    }

    template<class R, class X> inline
    void neg_(Interval<R>& r, const X& x) {
      neg_(r._lower,x);
      neg_(r._upper,x);
    }



    template<class R, class X, class Y> inline
    void add_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
      add_(r._lower,x._lower,y._lower,round_down); 
      add_(r._upper,x._upper,y._upper,round_up);
    }

    template<class R, class X, class Y> inline
    void add_(Interval<R>& r, const Interval<X>& x, const Y& y) {
      add_(r._lower,x._lower,y,round_down); 
      add_(r._upper,x._upper,y,round_up);
    }

    template<class R, class X, class Y> inline
    void add_(Interval<R>& r, const X& x, const Interval<Y>& y) {
      add_(r._lower,x,y._lower,round_down); 
      add_(r._upper,x,y._upper,round_up);
    }

    template<class R, class X, class Y> inline
    void add_(Interval<R>& r, const X& x, const Y& y) {
      add_(r._lower,x,y,round_down); 
      add_(r._upper,x,y,round_up);
    }



    template<class R, class X, class Y> inline
    void sub_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
      sub_(r._lower,x._lower,y._upper,round_down); 
      sub_(r._upper,x._upper,y._lower,round_up);
    }

    template<class R, class X, class Y> inline
    void sub_(Interval<R>& r, const Interval<X>& x, const Y& y) {
      sub_(r._lower,x._lower,y,round_down); 
      sub_(r._upper,x._upper,y,round_up);
    }

    template<class R, class X, class Y> inline
    void sub_(Interval<R>& r, const X& x, const Interval<Y>&  y) {
      sub_(r._lower,x,y._upper,round_down); 
      sub_(r._upper,x,y._lower,round_up);
    }

    template<class R, class X, class Y> inline
    void sub_(Interval<R>& r, const X& x, const Y& y) {
      sub_(r._lower,x,y,round_down); 
      sub_(r._upper,x,y,round_up);
    }




    

    



    template<class R, class X, class Y> inline
    void mul_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
      R& rl = r._lower;
      R& ru = r._upper;
      const X& xl = x.lower();
      const X& xu = x.upper();
      const Y& yl = y.lower();
      const Y& yu = y.upper();
      RoundDown d;
      RoundUp u;

      if (xl>=0) {
        if (yl>=0) {
          mul_(rl,xl,yl,d); mul_(ru,xu,yu,u);
        } else if(yu<=0) {
          mul_(rl,xu,yl,d); mul_(ru,xl,yu,u);
        } else {
          mul_(rl,xu,yl,d); mul_(ru,xu,yu,u);
        }
      } else if (xu<=0) {
        if (yl>=0) {
          mul_(rl,xl,yu,d); mul_(ru,xu,yl,u);
        } else if(yu<=0) {
          mul_(rl,xu,yu,d); mul_(ru,xl,yl,u); 
        } else {
          mul_(rl,xl,yu,d); mul_(ru,xl,yl,u);
        }
      } else {
        if (yl>=0) {
          mul_(rl,xl,yu,d); mul_(ru,xu,yu,u);
        } else if(yu<=0) {
          mul_(rl,xu,yl,d); mul_(ru,xl,yl,u);
        } else {
          R t1; R t2; 
          mul_(t1,xl,yu,d); mul_(t2,xu,yl,d); min_(rl,t1,t2);
          mul_(t1,xl,yl,u); mul_(t2,xu,yu,u); max_(ru,t1,t2);
        }
      }
    }    

    template<class R, class X, class Y> inline
    void mul_(Interval<R>& r, const Interval<X>& x, const Y& y) {
      mul_(r,x,Interval<R>(y));
    }

    template<class R, class X, class Y> inline
    void mul_(Interval<R>& r, const X& x, const Interval<Y>& y) {
      mul_(r,Interval<R>(x),y);
    }

    template<class R, class X, class Y> inline
    void mul_(Interval<R>& r, const X& x, const Y& y) {
      mul_(r._lower,x,y,round_down); mul_(r._upper,x,y,round_up);
    }


    template<class R, class X, class Y> inline
    void div_(Interval<R>& r, const Interval<X>& x, const Interval<Y>& y) {
      typedef Interval<R> I;
      R& rl = r._lower;
      R& ru = r._upper;
      const X& xl = x.lower();
      const X& xu = x.upper();
      const Y& yl = y.lower();
      const Y& yu = y.upper();
      RoundDown d;
      RoundUp u;

      if (yl>0) {
        if (xl>=0) {
          div_(rl,xl,yu,d); div_(ru,xu,yl,u);
        } else if(xu<=0) {
          div_(rl,xl,yl,d); div_(ru,xu,yu,u);
        } else {
          div_(rl,xl,yl,d); div_(ru,xu,yl,u);
        }
      } else if (yu<0) {
        if (xl>=0) {
          div_(rl,xu,yu,d); div_(ru,xl,yl,u);
        } else if(xu<=0) {
          div_(rl,xu,yl,d); div_(ru,xl,yu,u);
        } else {
          div_(rl,xu,yu,d); div_(ru,xl,yu,u);
        }
      } else {
        inf_(ru); neg_(rl,ru);
      }
    }

    template<class R, class X, class Y> inline
    void div_(Interval<R>& r, const Interval<X>& x, const Y& y) {
      div_(r,x,Interval<R>(y)); 
    }

    template<class R, class X, class Y> inline
    void div_(Interval<R>& r, const X& x, const Interval<Y>& y) {
      div_(r,Interval<R>(x),y);
    }

    template<class R, class X, class Y> inline
    void div_(Interval<R>& r, const X& x, const Y& y) {
      div_(r._lower,x,y,round_down);
      div_(r._upper,x,y,round_up);
    }




    template<class R> inline
    void min_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
      r._lower=min(x1.lower(),x2.lower()); r._upper=min(x1.upper(),x2.upper());
    }
    
    template<class R> inline
    void max_(Interval<R>& r, const Interval<R>& x1, const Interval<R>& x2) {
      r._lower=max(x1.lower(),x2.lower()); r._upper=max(x1.upper(),x2.upper());
    }


    template<class R> inline
    void abs_(Interval<R>& r, const Interval<R>& x) {
      if(x.lower()>=0) { r=x; } 
      if(x.upper() < 0) { r=-x; } 
      R nxl; neg_(nxl,x.lower()); r=Interval<R>(R(0),max(nxl,x.upper()));
    }
  
    template<class R, class X> inline
    void abs_(Interval<R>& r, const Interval<X>& x) {
      if(x.lower()>=0) { r=x; } 
      if(x.upper() < 0) { r=-x; } 
      R nxl; neg_(nxl,x.lower(),round_up); 
      R xu; pos_(xu,x.upper(),round_up);
      set_(r._lower,0); max_(r._upper,nxl,xu);
    }
  
    template<class R, class X> inline
    void abs_(Interval<R>& r, const X& x) {
      if(x>=0) { r=x; } else { r=-x; }
    }
  
    template<class R> inline
    Interval<R> abs(const Interval<R>& x) {
      Interval<R> r; abs_(r,x); return r; }

    template<class R,class N> inline
    void pow_(Interval<R>& r, const Interval<R>& x, const N& n) {
      Interval<R> result=R(1);
      for(N i=0; i!=n; ++i) {
        result*=x;
      }
      r=result;
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
    bool overlap(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.upper()>x2.lower() && x1.lower()<x2.upper());
    }

    template<class R> inline
    bool subset(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.lower()>=x2.lower() && x1.upper()<=x2.upper());
    }

    template<class R> inline
    bool inside(const Interval<R>& x1, const Interval<R>& x2) {
      return (x1.lower()>x2.lower() && x1.upper()<x2.upper());
    }

    template<class R> inline
    Interval<R> intersection(const Interval<R>& x1, const Interval<R>& x2) {
      return Interval<R>(max(x1.lower(),x2.lower()),
                         min(x1.upper(),x2.upper()));
    }

    template<class R> inline
    Interval<R> hull(const Interval<R>& x, const Interval<R>& y) {
      if(x.empty()) { return y; }
      if(y.empty()) { return x; }
      return Interval<R>(min(x.lower(),y.lower()),
                         max(x.upper(),y.upper()));
    }
    
   

    
    


   
    template<class R> template<class E> inline
    Interval<R>::Interval(const Expression<E>& e) { e.assign_to(*this); }

    template<class R> template<class E> inline
    Interval<R>& Interval<R>::operator=(const Expression<E>& e) { 
      e.assign_to(*this); return *this; }




    template<class R> 
    std::ostream& operator<<(std::ostream& os, const Interval<R>& ivl) 
    {
      if(ivl.empty()) {
        return os << "[1:0]";
      }
      else {
        return os << "[" << ivl.lower() << ":" << ivl.upper() << "]";
      }
    }
    
    
    template<class R> 
    std::istream& 
    operator>>(std::istream& is, Interval<R>& ivl)
    {
      char c;
      R l;
      R u;
      is >> c;
      if(c=='[') {
        is >> l >> c;
        if(c!=',' && c!=':') {
          is.setstate(std::ios_base::failbit);
        }
        is >> u >> c;
        if(c!=']') {
          is.setstate(std::ios_base::failbit);
        }
        ivl=Interval<R>(l,u);
      } else {
        is.putback(c);
        is >> l;
        ivl=Interval<R>(l);
      }
      return is;
    }
        
    
    
    
  } // namespace Numeric
} // namespace Ariadne
  
#include "numeric/transcendental.h"

namespace TBLAS {

  template<class real>
  int
  iamax_ (const int N, const Ariadne::Numeric::Interval<real> *X, const int incX)
  {
#ifdef DEBUG
    std::cerr << "TBLAS::iamax_(const int N, const interval<real> *X, const int incX)\n";
#endif

    real mx = 0;
    int ix = 0;
    int i;
    int result = 0;
    
    if (incX <= 0) {
      return 0;
    }

    for (i = 0; i < N; i++) {
      real av=Ariadne::Numeric::min_(
                Ariadne::Numeric::abs_(X[ix].lower()),
                Ariadne::Numeric::abs_(X[ix].upper()));
      if (av > mx) {
        mx = av;
        result = i;
      }
      ix += incX;
    }
    
    return result;
  }
  
}
