/***************************************************************************
 *            interval.code.h
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


template<class R>
Interval<R>::Interval(const std::string& str) {
  std::stringstream ss(str); ss>>*this; }

template<class R> void instantiate_interval();

template<typename R> void cos_(Interval<R>& res, const Interval<R>& ivl);



template<class R> 
void 
pi_(Interval<R>& r) 
{
  pi_(r._lower,round_down);
  pi_(r._upper,round_up);
}

template<class R> 
void 
sin_(Interval<R>& r, const Interval<R>& x) 
{
  Interval<R> y=x-pi< Interval<R> >()/2;
  cos_(r,y);
}



template<class R> 
void tan_(Interval<R>& r, const Interval<R>& x) 
{
  tan_(r._lower,x._lower,round_down);
  tan_(r._upper,x._upper,round_up);
}


template<class R> 
void asin_(Interval<R>& r, const Interval<R>& x) 
{
  asin_(r._lower,x._lower,round_down);
  asin_(r._upper,x._upper,round_up);
}


template<class R> 
void acos_(Interval<R>& r, const Interval<R>& x) 
{
  // Need to consider case r=x;
  Interval<R> t(x);
  acos_(r._lower,t._upper,round_down);
  acos_(r._upper,t._lower,round_up);
}


template<class R> 
void atan_(Interval<R>& r, const Interval<R>& x) 
{
  atan_(r._lower,x._lower,round_down);
  atan_(r._upper,x._upper,round_up);
}




template<typename R>  
void
cos_(Interval<R>& res, const Interval<R>& ivl) 
{
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
    //if(!(l>=0 && l<=tpu)) {
    //  std::cerr << "ivl=" << ivl << ", [l:u]=[" << l << ":" << u << "]" << std::endl;
    //}
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
    
/*

template<typename R>  
Interval<R> 
tan(const Interval<R>& ivl) 
{
  return sin(ivl)/cos(ivl);
}
    

template<typename R>  
Interval<R> 
asin(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: asin(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(asin_down(ivl.lower()),asin_up(ivl.upper()));
}


template<typename R>
Interval<R> 
acos(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: acos(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(acos_down(ivl.upper()),acos_up(ivl.lower()));
}
   
 
template<typename R>  
Interval<R> 
atan(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: atan(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(atan_down(ivl.lower()),atan_up(ivl.upper()));
}
    
*/

template<typename R>  
void
instantiate_interval()
{
  Interval<R>* ivl=0;
  const Interval<R>* civl=0;
  sqrt(*ivl);
  exp(*ivl);
  log(*ivl);
  pi_(*ivl);
  sin_(*ivl,*civl);
  sin(*ivl);
  cos_(*ivl,*civl);
  cos(*ivl);
  tan(*ivl);
  asin(*ivl);
  acos_(*ivl,*civl);
  acos(*ivl);
  atan(*ivl);
}



} // namespace Numeric
} // namespace Ariadne
  
