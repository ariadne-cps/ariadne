/***************************************************************************
 *            affine_model.code.h
 *
 *  Copyright  2007  Pieter Collins
 *  Pieter.Collins@cwi.nl
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
 
#include "affine_model.h"

#include "differentiation/affine_variable.h"
#include "function/function_interface.h"

#include "output/logging.h"

namespace Ariadne {

template<class R> 
void 
AffineModel<R>::instantiate()
{
  AffineModel<R>* am=0;
  compose(*am,*am);
  inverse(*am);
  implicit(*am);
}

template<class R> 
AffineModel<R>::AffineModel(const Vector<I>& d,
                            const Vector<R>& c, 
                            const Vector<I>& v, 
                            const Matrix<I>& j)
  : _domain(d), _centre(c), _value(v), _jacobian(j)
{
  ARIADNE_ASSERT(c.size()==d.size());
  ARIADNE_ASSERT(v.size()==j.number_of_rows());
  ARIADNE_ASSERT(c.size()==j.number_of_columns());
}

template<class R> 
AffineModel<R>::AffineModel(const Vector<I>& d,
                            const Vector<R>& c, 
                            const Vector<AffineVariable<I> >& av)
  : _domain(d), _centre(c), _value(av.size()), _jacobian(av.size(),c.size())
{
  ARIADNE_ASSERT(c.size()==d.size());
  for(uint i=0; i!=av.size(); ++i) {
    _value[i]=av[i].value();
    for(uint j=0; j!=c.size(); ++j) {
      _jacobian(i,j)=av[i].gradient(j);
    }
  }
}



template<class R> 
AffineModel<R>::AffineModel(const Vector<I>& d,
                            const Vector<R>& c,
                            const FunctionInterface<R>& f)
  : _domain(d), _centre(c), _value(f.evaluate(c)), _jacobian(f.jacobian(d))
{
}


template<class R>
Vector<typename traits<R>::interval_type> 
AffineModel<R>::evaluate(const Vector<I>& pt) const
{
  return this->_value+this->_jacobian*(pt-this->_centre);
}



/*

template<class R> 
AffineModel<R>
operator+(const AffineModel<R>& am) 
{
  return am;
}

template<class R> 
AffineModel<R>
operator-(const AffineModel<R>& am) 
{
  return AffineModel<R>(am.domain(),am.centre(),-am.value(),-am.jacobian());
}


template<class R> 
AffineModel<R>
operator+(const AffineModel<R>& am1, const AffineModel<R>& am2) 
{
  ARIADNE_LOG(3,"operator+(AffineModel am1, AffineModel am2)\n");
  if(am1.domain()==am2.domain() && am1.centre()==am2.centre()) {
    return AffineModel<R>(am1.domain(),am1.centre(),am1.value()+am2.value(),am1.jacobian()+am2.jacobian());
  } else {
    Box<R> nd=open_intersection(am1.domain(),am2.domain());
    return restrict(am1,nd)+restrict(am1,nd);
  }
}

template<class R> 
AffineModel<R>
operator-(const AffineModel<R>& am1, const AffineModel<R>& am2) 
{
  ARIADNE_LOG(3,"operator+(AffineModel am1, AffineModel am2)\n");
  if(am1.domain()==am2.domain() && am1.centre()==am2.centre()) {
    return AffineModel<R>(am1.domain(),am1.centre(),am1.value()-am2.value(),am1.jacobian()-am2.jacobian());
  } else {
    Box<R> nd=open_intersection(am1.domain(),am2.domain());
    return restrict(am1,nd)-restrict(am1,nd);
  }
}

*/


template<class R> 
Vector<typename AffineModel<R>::I>
AffineModel<R>::range() const
{
  return this->evaluate(this->_domain);
}


template<class R> 
AffineModel<R>
translate(const AffineModel<R>& am, const Vector<R>& nc) 
{
  typedef Interval<R> I;
  Vector<I> nd=am.domain();
  Vector<I> nv=am.evaluate(nc);
  const Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
AffineModel<R>
restrict(const AffineModel<R>& am, const Vector< Interval<R> >& nd) 
{
  typedef Interval<R> I;
  Vector<R> nc=nd.centre();
  Vector<I> nv=am.evaluate(nc);
  const Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
AffineModel<R>
compose(const AffineModel<R>& am1, const AffineModel<R>& am2) 
{
  ARIADNE_ASSERT(refines(am2.range(),am1.domain()));
  typedef Interval<R> I;
  Vector<I> const& nd=am2.domain();
  Vector<R> const& nc=am2.centre();
  Vector<I> nv=am1.evaluate(am2.value());
  Matrix<I> nj=am1.jacobian()*am2.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}




template<class R> 
AffineModel<R>
reduce(const AffineModel<R>& am, size_type s) 
{
  assert(s<=am._s);
  if(s==am._s) { return am;  }
  
  const size_type rd=am._rd;
  const size_type ad=am._ad;
  
  AffineModel<R> res(rd,ad,s);
  typename AffineModel<R>::I tmp;
  for(size_type i=0; i!=rd; ++i) {
    res._iptr[i]=am._iptr[i];
    for(size_type j=i+rd; j!=i+rd*(ad+1); j+=rd) {
      res._rptr[j]=am._iptr[j].midpoint();
      res._iptr[i]+=(am._iptr[j]-res._rptr[j]);
    }
  }
  return res;
}




template<class R> 
AffineModel<R> 
inverse(const AffineModel<R>&) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> 
AffineModel<R> 
implicit(const AffineModel<R>&)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> 
std::ostream& 
AffineModel<R>::write(std::ostream& os) const
{
  os << "AffineModel"
     << "( domain=" << this->domain()
     << ", centre=" << this->centre()
     << ", value=" << this->value()
     << ", jacobian=" << this->jacobian()
     << " )";
  return os;
}



} // namespace Ariadne
