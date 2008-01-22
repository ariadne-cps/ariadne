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
Function::AffineModel<R>::instantiate()
{
  AffineModel<R>* am=0;
  compose(*am,*am);
  inverse(*am);
  implicit(*am);
}

template<class R> 
Function::AffineModel<R>::AffineModel(const Geometry::Box<R>& d,
                                      const Geometry::Point<R>& c, 
                                      const Geometry::Point<I>& v, 
                                      const LinearAlgebra::Matrix<I>& j)
  : _domain(d), _centre(c), _value(v), _jacobian(j)
{
  ARIADNE_ASSERT(c.dimension()==d.dimension());
  ARIADNE_ASSERT(v.dimension()==j.number_of_rows());
  ARIADNE_ASSERT(c.dimension()==j.number_of_columns());
}

template<class R> 
Function::AffineModel<R>::AffineModel(const Geometry::Box<R>& d,
                                      const Geometry::Point<R>& c, 
                                      const array<Function::AffineVariable<I> >& av)
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
Function::AffineModel<R>::AffineModel(const Geometry::Box<R>& d,
                                      const Geometry::Point<R>& c,
                                      const Function::FunctionInterface<R>& f)
  : _domain(d), _centre(c), _value(f.evaluate(c.position_vector())), _jacobian(f.jacobian(d.position_vectors()))
{
}


template<class R>
Geometry::Point<typename Numeric::traits<R>::interval_type> 
Function::AffineModel<R>::evaluate(const Geometry::Point<I>& pt) const
{
  return this->_value+this->_jacobian*(pt-this->_centre);
}



/*

template<class R> 
Function::AffineModel<R>
Function::operator+(const Function::AffineModel<R>& am) 
{
  return am;
}

template<class R> 
Function::AffineModel<R>
Function::operator-(const Function::AffineModel<R>& am) 
{
  return AffineModel<R>(am.domain(),am.centre(),-am.value(),-am.jacobian());
}


template<class R> 
Function::AffineModel<R>
Function::operator+(const Function::AffineModel<R>& am1, const Function::AffineModel<R>& am2) 
{
  ARIADNE_LOG(3,"operator+(Function::AffineModel am1, Function::AffineModel am2)\n");
  if(am1.domain()==am2.domain() && am1.centre()==am2.centre()) {
    return AffineModel<R>(am1.domain(),am1.centre(),am1.value()+am2.value(),am1.jacobian()+am2.jacobian());
  } else {
    Geometry::Box<R> nd=open_intersection(am1.domain(),am2.domain());
    return restrict(am1,nd)+restrict(am1,nd);
  }
}

template<class R> 
Function::AffineModel<R>
Function::operator-(const Function::AffineModel<R>& am1, const Function::AffineModel<R>& am2) 
{
  ARIADNE_LOG(3,"operator+(Function::AffineModel am1, Function::AffineModel am2)\n");
  if(am1.domain()==am2.domain() && am1.centre()==am2.centre()) {
    return AffineModel<R>(am1.domain(),am1.centre(),am1.value()-am2.value(),am1.jacobian()-am2.jacobian());
  } else {
    Geometry::Box<R> nd=open_intersection(am1.domain(),am2.domain());
    return restrict(am1,nd)-restrict(am1,nd);
  }
}

*/


template<class R> 
Function::AffineModel<R>
Function::translate(const Function::AffineModel<R>& am, const Geometry::Point<R>& nc) 
{
  typedef Numeric::Interval<R> I;
  Geometry::Box<R> nd=am.domain();
  Geometry::Point<I> nv=am.evaluate(nc);
  const LinearAlgebra::Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
Function::AffineModel<R>
Function::restrict(const Function::AffineModel<R>& am, const Geometry::Box<R>& nd) 
{
  typedef Numeric::Interval<R> I;
  Geometry::Point<R> nc=nd.centre();
  Geometry::Point<I> nv=am.evaluate(nc);
  const LinearAlgebra::Matrix<I> & nj=am.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}

template<class R> 
Function::AffineModel<R>
Function::compose(const Function::AffineModel<R>& am1, const Function::AffineModel<R>& am2) 
{
  ARIADNE_ASSERT(Geometry::subset(am2.range(),am1.domain()));
  typedef Numeric::Interval<R> I;
  Geometry::Box<R> const& nd=am2.domain();
  Geometry::Point<R> const& nc=am2.centre();
  Geometry::Point<I> nv=am1.evaluate(reinterpret_cast<const Geometry::Point<I>&>(am2.value()));
  LinearAlgebra::Matrix<I> nj=am1.jacobian()*am2.jacobian();
  return AffineModel<R>(nd,nc,nv,nj);
}




template<class R> 
Function::AffineModel<R>
Function::reduce(const Function::AffineModel<R>& am, size_type s) 
{
  assert(s<=am._s);
  if(s==am._s) { return am;  }
  
  const dimension_type rd=am._rd;
  const dimension_type ad=am._ad;
  
  Function::AffineModel<R> res(rd,ad,s);
  typename Function::AffineModel<R>::I tmp;
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
Function::AffineModel<R> 
Function::inverse(const Function::AffineModel<R>&) 
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> 
Function::AffineModel<R> 
Function::implicit(const Function::AffineModel<R>&)
{
  throw NotImplemented(__PRETTY_FUNCTION__);
}

template<class R> 
std::ostream& 
Function::AffineModel<R>::write(std::ostream& os) const
{
  os << "AffineModel"
     << "( domain=" << this->domain()
     << ", centre=" << this->centre()
     << ", value=" << this->value()
     << ", jacobian=" << this->jacobian()
     << " )";
  return os;
}






}
