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
 
namespace Ariadne {

namespace Numeric { template<class R> void instantiate_interval(); }

template<typename R>  
Numeric::Interval<R> 
Numeric::sin(const Interval<R>& ivl) 
{
  if(ivl.lower()>-1.5 && ivl.upper()<+1.5) {
    return Interval<R>(sin_down(ivl.lower()),sin_up(ivl.upper()));
  } else {
    return cos(ivl-pi<R>()/2);
  }
}
    

template<typename R>  
Numeric::Interval<R> 
Numeric::cos(const Interval<R>& ivl) 
{
  assert(ivl.lower()<=ivl.upper());
  Interval<R> pi=Numeric::pi<R>();
  Interval<R> two_pi=2*pi;
  R n=floor(div_down(ivl.lower(),two_pi.lower()));
  //std::cerr << std::setprecision(20);
  //std::cerr << "n=" << n << std::endl;       
  // l and u are intervals shifted by a multiple of 2*pi.
  Interval<R> x=ivl-n*two_pi;
  const R& pl=pi.lower();
  const R& pu=pi.upper();
  const R& tpl=two_pi.lower();
  const R& tpu=two_pi.upper();
  const R  thpl=mul_down(pi.lower(),3);
  const R& l=x.lower();
  const R& u=x.upper();
  //std::cerr << "l=" << l << ", u=" << u << std::endl;       
  if(sub_down(u,l)>=tpu) {
    // There is a full period in [l,u]
    return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
  } else {
    // Check values are reasonable
    //if(!(l>=0 && l<=tpu)) {
    //  std::cerr << "ivl=" << ivl << ", [l:u]=[" << l << ":" << u << "]" << std::endl;
    //}
    assert(l>=-0 && l<=tpu);
    if(l<pl) {
      if(u<pl) {
        return Interval<R>(cos_down(u),cos_up(l));
      } else if(u<tpl) {
        return Interval<R>(static_cast<R>(-1),max(cos_up(l),cos_up(u)));
      } else {
        return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
      }
    } else if (l<pu) {
      if(u<tpl) {
        return Interval<R>(static_cast<R>(-1),max(cos_up(l),cos_up(u)));
      } else {
        return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
      }
    } else if (l<tpl) {
      if(u<tpl) {
        return Interval<R>(cos_down(l),cos_up(u));
      } else if(u<thpl) {
        return Interval<R>(min(cos_down(l),cos_down(u)),static_cast<R>(1));
      } else {
        return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
      }
    } else { // 2*pi_ < ln < 2*pi^
      if(u<thpl) {
            return Interval<R>(min(cos_down(l),cos_down(u)),static_cast<R>(1));
      } else {
        return Interval<R>(static_cast<R>(-1),static_cast<R>(1));
      }
    }
  }
}
    

template<typename R>  
Numeric::Interval<R> 
Numeric::tan(const Interval<R>& ivl) 
{
  return sin(ivl)/cos(ivl);
}
    

template<typename R>  
Numeric::Interval<R> 
Numeric::asin(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: asin(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(asin_down(ivl.lower()),asin_up(ivl.upper()));
}


template<typename R>
Numeric::Interval<R> 
Numeric::acos(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: acos(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(acos_down(ivl.upper()),acos_up(ivl.lower()));
}
   
 
template<typename R>  
Numeric::Interval<R> 
Numeric::atan(const Interval<R>& ivl) 
{
  static bool warn=true;
  if(warn) {
    std::cerr << "WARNING: atan(Interval) does not handle branch cuts correctly" << std::endl;
    warn=false;
  }
  return Interval<R>(atan_down(ivl.lower()),atan_up(ivl.upper()));
}
    
template<typename R>  
void
Numeric::instantiate_interval()
{
  Interval<R>* ivl=0;
  sin(*ivl);
  cos(*ivl);
  tan(*ivl);
  asin(*ivl);
  acos(*ivl);
  atan(*ivl);
}

} // namespace Ariadne
  
