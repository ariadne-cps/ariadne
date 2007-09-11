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
 
namespace Ariadne {
  namespace Numeric {
  
    template<class R> inline 
    Interval<R>::Interval()
      : _lower(conv_down<R>(0)), _upper(conv_up<R>(0)) { }

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
      : _lower(conv_down<R>(l)), _upper(conv_up<R>(u)) { }

    template<class R> template<class RX> inline 
    Interval<R>::Interval(const Interval<RX>& ivl)
      : _lower(conv_down<R>(ivl.lower())), _upper(conv_up<R>(ivl.upper())) { }

    template<class R> template<class RX> inline 
    Interval<R>::Interval(const RX& x)
      : _lower(conv_down<R>(x)), _upper(conv_up<R>(x)) { }
      
    template<class R> template<class RX> inline 
    Interval<R>& Interval<R>::operator=(const RX& x) {
      this->_lower=conv_down<R>(x); this->_upper=conv_up<R>(x); return *this;
    }
  

    template<class R> template<class RX> inline 
    Interval<R>& Interval<R>::operator=(const Interval<RX>& ivl) {
      this->_lower=conv_down<R>(ivl.lower()); this->_upper=conv_up<R>(ivl.upper()); return *this;
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
      return div_approx(add_approx(this->lower(),this->upper()),R(2));
    }

    template<class R> inline 
    R Interval<R>::radius() const { 
      return div_up(sub_up(this->upper(),this->lower()),R(2)); 
    }

    template<class R> inline
    R Interval<R>::width() const { 
      return sub_up(this->upper(),this->lower()); 
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
    void Interval<R>::expand_by(const R& r) { 
      this->_lower=sub_down(this->lower(),r); this->_upper=add_up(this->upper(),r); 
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


    template<class R, class RX> inline
    bool encloses(const Interval<R>& ivl, const RX& x) { 
      return ivl.encloses(x);
    }

    template<class R1, class R2> inline
    bool refines(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      return ivl1.refines(ivl2);
    }


    template<class R1, class R2> inline
    tribool operator==(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()==ivl2.upper() && ivl1.upper()==ivl2.lower()) { return true; }
      if(ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower()) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator!=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>ivl2.upper() || ivl1.upper()<ivl2.lower()) { return true; }
      if(ivl1.lower()==ivl2.upper() && ivl1.upper()==ivl2.lower()) { return false; }
      return indeterminate;
    }
        
    
    template<class R1, class R2> inline
    tribool operator<(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.upper()<ivl2.lower()) { return true; } 
      if(ivl1.lower()>=ivl2.upper()) { return false; }
      return indeterminate;
    }
        
    
    template<class R1, class R2> inline
    tribool operator>(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>ivl2.upper()) { return true; }
      if(ivl1.upper()<=ivl2.lower()) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator<=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.upper()<=ivl2.lower()) { return true; } 
      if(ivl1.lower()>ivl2.upper()) { return false; }
      return indeterminate;
    }
        
    
    template<class R1, class R2> inline
    tribool operator>=(const Interval<R1>& ivl1, const Interval<R2>& ivl2) { 
      if(ivl1.lower()>=ivl2.upper()) { return true; }
      if(ivl1.upper()<ivl2.lower()) { return false; }
      return indeterminate;
    }
      


    
    template<class R1, class R2> inline
    tribool operator==(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()==x && ivl.upper()==x) { return true; }
      if(ivl.lower()>x || ivl.upper()<x) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator!=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>x || ivl.upper()<x) { return true; }
      if(ivl.lower()==x && ivl.upper()==x) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator<(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.upper()<x) { return true; }
      if(ivl.lower()>=x) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator>(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>x) { return true; }
      if(ivl.upper()<=x) { return false; }
      return indeterminate;
    }
      
    
    template<class R1, class R2> inline
    tribool operator<=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.upper()<=x) { return true; }
      if(ivl.lower()>x) { return false; }
      return indeterminate;
    }

    
    template<class R1, class R2> inline
    tribool operator>=(const Interval<R1>& ivl, const R2& x) { 
      if(ivl.lower()>=x) { return true; }
      if(ivl.upper()<x) { return false; }
      return indeterminate;
    }

    
    template<class R1, class R2> inline
    tribool operator==(const R1& x, const Interval<R2>& ivl) { 
      return ivl==x; 
    }

      
    
    template<class R1, class R2> inline
    tribool operator!=(const R1& x, const Interval<R2>& ivl) {
      return ivl!=x;
    }
      
    
    template<class R1, class R2> inline
    tribool operator<(const R1& x, const Interval<R2>& ivl) {
      return ivl>x;
    }
      
    
    template<class R1, class R2> inline
    tribool operator>(const R1& x, const Interval<R2>& ivl) {
      return ivl<x;
    }

    
    template<class R1, class R2> inline
    tribool operator<=(const R1& x, const Interval<R2>& ivl) {
      return ivl>=x;
    }
      
    
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
    Interval<R> operator+(const Interval<R>& x1, const int& x2) {
      return x1+Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator-(const Interval<R>& x1, const int& x2) {
      return x1-Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator*(const Interval<R>& x1, const int& x2) {
      return x1*Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator/(const Interval<R>& x1, const int& x2) {
      return x1/Interval<R>(x2); 
    }


    template<class R> inline
    Interval<R> operator+(const int& x1, const Interval<R>& x2) {
      return Interval<R>(x1)+x2;
    }

    template<class R> inline
    Interval<R> operator-(const int& x1, const Interval<R>& x2) {
      return Interval<R>(x1)-Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator*(const int& x1, const Interval<R>& x2) {
      return Interval<R>(x1)*Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator/(const int& x1, const Interval<R>& x2) {
      return Interval<R>(x1)/Interval<R>(x2); 
    }

    
    
    template<class R> inline
    Interval<R> operator+(const Interval<R>& x1, const double& x2) {
      return x1+Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator-(const Interval<R>& x1, const double& x2) {
      return x1-Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator*(const Interval<R>& x1, const double& x2) {
      return x1*Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator/(const Interval<R>& x1, const double& x2) {
      return x1/Interval<R>(x2); 
    }

    template<class R> inline
    Interval<R> operator+(const double& x1, const Interval<R>& x2) {
      return Interval<R>(x1)+x2; 
    }

    template<class R> inline
    Interval<R> operator-(const double& x1, const Interval<R>& x2) {
      return Interval<R>(x1)-x2;
    }

    template<class R> inline
    Interval<R> operator*(const double& x1, const Interval<R>& x2) {
      return Interval<R>(x1)*x2;
    }

    template<class R> inline
    Interval<R> operator/(const double& x1, const Interval<R>& x2) {
      return Interval<R>(x1)/x2; 
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
    Interval<R> pi() {
      return Interval<R>(pi_down<R>(),pi_up<R>());
    }


    template<typename R> Interval<R> sin(const Interval<R>& ivl);
    template<typename R> Interval<R> cos(const Interval<R>& ivl);
    template<typename R> Interval<R> tan(const Interval<R>& ivl);
    template<typename R> Interval<R> asin(const Interval<R>& ivl);
    template<typename R> Interval<R> acos(const Interval<R>& ivl);
    template<typename R> Interval<R> atan(const Interval<R>& ivl);
    
    
    template<class R> 
    std::ostream& Interval<R>::write(std::ostream& os) const {
      if(this->empty()) {
        return os << "[1:0]";
      }
      else {
        return os << "[" << this->lower() << ":" << this->upper() << "]";
      }
    }
    
    
    template<class R> 
    std::istream& Interval<R>::read(std::istream& is) {
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
        (*this)=Interval<R>(l,u);
      } else {
        is.putback(c);
        is >> l;
        (*this)=Interval<R>(l);
      }
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
