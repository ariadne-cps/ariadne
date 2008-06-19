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
 
#include "affine_variable.h"
#include "affine_model.h"

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
AffineModel<R>::AffineModel(const Box<R>& d,
                                      const Point<R>& c, 
                                      const Point<I>& v, 
                                      const Matrix<I>& j)
  : _domain(d), _centre(c), _value(v), _jacobian(j)
{
  ARIADNE_ASSERT(c.dimension()==d.dimension());
  ARIADNE_ASSERT(v.dimension()==j.number_of_rows());
  ARIADNE_ASSERT(c.dimension()==j.number_of_columns());
}

template<class R> 
AffineModel<R>::AffineModel(const Box<R>& d,
                                      const Point<R>& c, 
                                      const array<AffineVariable<I> >& av)
  : _domain(d), _centre(c), _value(av.size()), _jacobian(av.size(),c.dimension())
{
  ARIADNE_ASSERT(c.dimension()==d.dimension());
  for(uint i=0; i!=av.size(); ++i) {
    _value[i]=av[i].value();
    for(uint j=0; j!=c.dimension(); ++j) {
      _jacobian(i,j)=av[i].gradient(j);
    }
  }
}



template<class R> 
AffineModel<R>::AffineModel(const Box<R>& d,
                                      const Point<R>& c,
                                      const FunctionInterface<R>& f)
  : _domain(d), _centre(c), _value(f.evaluate(c.position_vector())), _jacobian(f.jacobian(d.position_vectors()))
{
}


template<class R>
Point<typename traits<R>::interval_type> 
AffineModel<R>::evaluate(const Point<I>& pt) const
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
Box<R>
AffineModel<R>::range() const
{
  return Box<R>(this->evaluate(Point<I>(this->_domain.position_vectors())));
}


template<class R> 
AffineModel<R>
translate(const AffineModel<R>& am, const Point<R>& nc) 
{
  typedef Interval<R> I;
  Box<R> nd=am.domain();
  Point<I> nv=am.evaluate(nc);
  const Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
AffineModel<R>
restrict(const AffineModel<R>& am, const Box<R>& nd) 
{
  typedef Interval<R> I;
  Point<R> nc=nd.centre();
  Point<I> nv=am.evaluate(nc);
  const Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
AffineModel<R>
compose(const AffineModel<R>& am1, const AffineModel<R>& am2) 
{
  ARIADNE_ASSERT(subset(am2.range(),am1.domain()));
  typedef Interval<R> I;
  Box<R> const& nd=am2.domain();
  Point<R> const& nc=am2.centre();
  Point<I> nv=am1.evaluate(reinterpret_cast<const Point<I>&>(am2.value()));
  Matrix<I> nj=am1.jacobian()*am2.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}




template<class R> 
AffineModel<R>
reduce(const AffineModel<R>& am, size_type s) 
{
  assert(s<=am._s);
  if(s==am._s) { return am;  }
  
  const dimension_type rd=am._rd;
  const dimension_type ad=am._ad;
  
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
