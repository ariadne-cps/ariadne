/***************************************************************************
 *            interval-float.code.h
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
 
#include <sstream>
#include "numeric/interval.h"

namespace Ariadne {
namespace Numeric { 


template<class R> void instantiate_interval();

template<typename T> void cos_(Interval< Float<T> >& res, const Interval< Float<T> >& ivl);



template<class T> 
void 
pi_(Interval< Float<T> >& r) 
{
  pi_(r._lower,round_down);
  pi_(r._upper,round_up);
}

template<class T> 
void 
sin_(Interval< Float<T> >& r, const Interval< Float<T> >& x) 
{
  Interval< Float<T> > y=x-pi< Interval< Float<T> > >()/2;
  cos_(r,y);
}



template<class T> 
void tan_(Interval< Float<T> >& r, const Interval< Float<T> >& x) 
{
  tan_(r._lower,x._lower,round_down);
  tan_(r._upper,x._upper,round_up);
}


template<class T> 
void asin_(Interval< Float<T> >& r, const Interval< Float<T> >& x) 
{
  asin_(r._lower,x._lower,round_down);
  asin_(r._upper,x._upper,round_up);
}


template<class T> 
void acos_(Interval< Float<T> >& r, const Interval< Float<T> >& x) 
{
  // Need to consider case r=x;
  Interval< Float<T> > t(x);
  acos_(r._lower,t._upper,round_down);
  acos_(r._upper,t._lower,round_up);
}


template<class T> 
void atan_(Interval< Float<T> >& r, const Interval< Float<T> >& x) 
{
  atan_(r._lower,x._lower,round_down);
  atan_(r._upper,x._upper,round_up);
}




template<class T>  
void
cos_(Interval< Float<T> >& res, const Interval< Float<T> >& ivl) 
{
  typedef Float<T> R;
  assert(ivl.lower()<=ivl.upper());
  Interval<R> pi=Numeric::pi< Interval<R> >();
  Interval<R> two_pi=2*pi;
  R n=floor<R>(div_down(ivl.lower(),two_pi.lower()));
  //std::cerr << std::setprecision(20);
  //std::cerr << "n=" << n << std::endl;       
  // l and u are intervals shifted by a mul_tiple of 2*pi.
  Interval<R> x=ivl-n*two_pi;
  const R& pl=pi.lower();
  const R& pu=pi.upper();
  const R& tpl=two_pi.lower();
  const R& tpu=two_pi.upper();
  const R  thpl=mul_down(pi.lower(),3);
  const R& l=x.lower();
  const R& u=x.upper();
  Interval<R>& r=res;
  //R& rl=r._lower;
  //R& ru=r._upper;
  //std::cerr << "l=" << l << ", u=" << u << std::endl;       
  if(sub_down(u,l)>=tpu) {
    // There is a full period in [l,u]
    set_(r._lower,-1);
    set_(r._upper,+1);
    return;
  } else {
    // Check values are reasonable
    if(!(l>=0 && l<=tpu)) {
      std::cerr << "ivl=" << ivl << ", [l:u]=[" << l << ":" << u << "]" << std::endl;
    }
    assert(l>=-0 && l<=tpu);
    if(l<pl) {
      if(u<pl) { 
        //Interval<R>(cos_down(u),cos_up(l));
        cos_(r._lower,x._lower,round_up);
        cos_(r._upper,x._upper,round_down);
        std::swap(r._lower,r._upper);
      } else if(u<tpl) {
        //Interval<R>(-1,max(cos_up(l),cos_up(u)));
        cos_(r._lower,x._lower,round_up);
        cos_(r._upper,x._upper,round_up);
        max_(r._upper,r._lower,r._upper);
        set_(r._lower,-1);
      } else {
        //Interval<R>(-1,1);
        set_(r._lower,-1);
        set_(r._upper,+1);
      }
    } else if (l<pu) {
      if(u<tpl) {
        //Interval<R>(static_cast<R>(-1),max(cos_up(l),cos_up(u)));
        cos_(r._upper,x._upper,round_up);
        cos_(r._lower,x._lower,round_up);
        max_(r._upper,r._lower,r._upper);
        set_(r._lower,-1);
      } else {
        //Interval<R>(static_cast<R>(-1),static_cast<R>(1));
        set_(r._lower,-1);
        set_(r._upper,+1);
      }
    } else if (l<tpl) {
      if(u<tpl) {
        //Interval<R>(cos_down(l),cos_up(u));
        cos_(r._lower,x._lower,round_down);
        cos_(r._upper,x._upper,round_up);
      } else if(u<thpl) {
        //Interval<R>(min(cos_down(l),cos_down(u)),+1);
        cos_(r._lower,x._lower,round_down);
        cos_(r._upper,x._upper,round_down);
        min_(r._lower,r._lower,r._upper);
        set_(r._upper,+1);
      } else {
        //Interval<R>(-1,+1);
        set_(r._lower,-1);
        set_(r._upper,+1);
      }
    } else { // 2*pi_ < ln < 2*pi^
      if(u<thpl) {
        //Interval<R>(min(cos_down(l),cos_down(u)),+1);
        cos_(r._lower,x._lower,round_down);
        cos_(r._upper,x._upper,round_down);
        min_(r._lower,r._lower,r._upper);
        set_(r._upper,+1);
      } else {
        //Interval<R>(-1,+1);
        set_(r._lower,-1);
        set_(r._upper,+1);
      }
    }
  }
}
    

template<class R>  
void
instantiate_interval()
{
  Interval<R>* ivl=0;
  const Interval<R>* civl=0;
  sqrt_(*ivl,*civl);
  exp_(*ivl,*civl);
  log_(*ivl,*civl);
  pi_(*ivl);
  sin_(*ivl,*civl);
  cos_(*ivl,*civl);
  tan_(*ivl,*civl);
  asin_(*ivl,*civl);
  acos_(*ivl,*civl);
  atan_(*ivl,*civl);
}



} // namespace Numeric
} // namespace Ariadne
  
