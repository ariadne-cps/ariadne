/***************************************************************************
 *            python/export_affine_function.cc
 *
 *  Copyright  2007  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, Pieter.Collins@cwi.nl
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

#include "python/float.h"

#include "numeric/rational.h"
#include "linear_algebra/vector.h"
#include "linear_algebra/covector.h"
#include "linear_algebra/matrix.h"
#include "differentiation/taylor_derivative.h"
#include "differentiation/sparse_differential.h"
#include "function/affine_function.h"
#include "function/identity_function.h"

#include "geometry/polyhedron.h"
#include "geometry/polytope.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;
using namespace Ariadne::Python;


//FIXME: This code should go elsewhere
template<class X>
Polyhedron<typename traits<X>::arithmetic_type> 
preimage(const AffineFunction<X>& af, const Polyhedron<X>& plhd)
{
  // Function f(x) = Ax+b = y
  const Matrix<X>& fA=af.A();
  const Vector<X>& fb=af.b();

  // Polyhedron Ay<=b
  const Matrix<X>& pA=plhd.A();
  const Vector<X>& pb=plhd.b();

  // Preimage is polyhedron pA fA x + pA fb <= pb
  return Polyhedron<typename traits<X>::arithmetic_type>(pA*fA,pb-pA*fb);
}


template<class X>
Polytope<typename traits<X>::arithmetic_type> 
image(const AffineFunction<X>& af, const Polytope<X>& pltp)\
{
  // Function f(x) = Ax+b
  const Matrix<X>& fA=af.A();
  const Vector<X>& fb=af.b();

  // Construct the transformation on rays
  Matrix<X> fT(fA.number_of_rows()+1, fA.number_of_columns()+1);
  fT(slice(0,fA.number_of_rows()),slice(0,fA.number_of_columns()))=fA;
  for(size_type i=0; i!=fb.size(); ++i) { fT(i,fA.number_of_columns())=fb[i]; }
  fT(fA.number_of_rows(),fA.number_of_columns())=1;
     
  // Polytope Vs
  const Matrix<X>& pV=pltp.generators();              

  // Preimage is polyhedron AV+o^T 
  return Polytope<X>(fT*pV);
}

template<class R>
void export_affine_function() 
{
  typedef typename traits<R>::arithmetic_type I;
  typedef typename traits<R>::approximate_arithmetic_type A;

  class_< AffineFunction<R>, bases< FunctionInterface<R> > > affine_function_class("AffineFunction",init< Matrix<I>,Vector<I> >());
  affine_function_class.def(init< Matrix<R>,Vector<R> >());
  affine_function_class.def(init< Vector<I>,Matrix<I> >());
  affine_function_class.def(init< Vector<R>,Matrix<R> >());
  affine_function_class.def("argument_size", &AffineFunction<R>::argument_size);
  affine_function_class.def("result_size", &AffineFunction<R>::result_size);
  affine_function_class.def("smoothness", &AffineFunction<R>::smoothness);
  affine_function_class.def("__call__",(Vector<I>(AffineFunction<R>::*)(const Vector<I>&)const)(&AffineFunction<R>::evaluate));
  affine_function_class.def("evaluate",(Vector<I>(AffineFunction<R>::*)(const Vector<I>&)const)(&AffineFunction<R>::evaluate));
  affine_function_class.def("jacobian",(Matrix<I>(AffineFunction<R>::*)(const Vector<I>&)const)(&AffineFunction<R>::jacobian));
  affine_function_class.def("derivative",(TaylorDerivative<I>(AffineFunction<R>::*)(const Vector<I>&)const)(&AffineFunction<R>::derivative));
  affine_function_class.def("expansion",(SparseDifferentialVector<A>(AffineFunction<R>::*)(const Vector<A>&, const ushort&)const)(&AffineFunction<R>::expansion));
  affine_function_class.def(self_ns::str(self));

  class_< IdentityFunction<R>, bases< FunctionInterface<R> > > identity_function_class("IdentityFunction",init<uint>());
  identity_function_class.def(self_ns::str(self));
}

template<>
void export_affine_function<Rational>() 
{
  typedef Rational Q;

  class_< AffineFunction<Q>, bases< FunctionInterface<Q> > > affine_function_class("QAffineFunction",init< Matrix<Q>,Vector<Q> >());
  affine_function_class.def(init< Matrix<Q>,Vector<Q> >());
  affine_function_class.def("argument_size", &AffineFunction<Q>::argument_size);
  affine_function_class.def("result_size", &AffineFunction<Q>::result_size);
  affine_function_class.def("smoothness", &AffineFunction<Q>::smoothness);
  affine_function_class.def("evaluate",(Vector<Q>(AffineFunction<Q>::*)(const Vector<Q>&)const)(&AffineFunction<Q>::evaluate));
  affine_function_class.def("jacobian",(Matrix<Q>(AffineFunction<Q>::*)(const Vector<Q>&)const)(&AffineFunction<Q>::jacobian));
  affine_function_class.def(self_ns::str(self));

  def("image",(Polytope<Q>(*)(const AffineFunction<Q>&,const Polytope<Q>&))&image);
  def("preimage",(Polyhedron<Q>(*)(const AffineFunction<Q>&,const Polyhedron<Q>&))&preimage);
}

template void export_affine_function<Rational>();
template void export_affine_function<FloatPy>();
