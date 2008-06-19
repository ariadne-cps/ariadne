/***************************************************************************
 *            standard_integrator.code.h
 *
 *  Copyright  2006-7  Alberto Casagrande, Pieter Collins
 *  casagrande@dimi.uniud.it, pieter.collins@cwi.nl
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
 
//#define DEBUG

#include <iosfwd>
#include <string>
#include <sstream>
#include <algorithm>
#include <typeinfo>

#include <list>
#include <set>
#include <vector>
#include <valarray>

#include "standard_integrator.h"

#include "base/array.h"
#include "base/exceptions.h"
#include "numeric/interval.h"

#include "linear_algebra/vector.h"

#include "linear_algebra/matrix.h"
#include "linear_algebra/matrix_function.h"

#include "geometry/box.h"
#include "geometry/zonotope.h"

#include "function/affine_model.h"
#include "function/taylor_model.h"

#include "system/vector_field.h"

#include "evaluation/standard_bounder.h"
#include "evaluation/standard_flower.h"

#include "output/logging.h"



namespace Ariadne { 


template<class R> 
StandardIntegrator< Zonotope<R> >::
StandardIntegrator() 
  : _bounder(new StandardBounder<R>()),
    _flower(new StandardFlower<R>(4u,1u))
{ }

template<class R> 
StandardIntegrator< Zonotope<R> >::
StandardIntegrator(const BounderInterface<R>& bounder, 
                   const FlowerInterface<R>& flower)
  : _bounder(bounder.clone()),
    _flower(flower.clone())
{ }

template<class R> 
StandardIntegrator< Zonotope<R> >*
StandardIntegrator< Zonotope<R> >::
clone() const
{
  return new StandardIntegrator< Zonotope<R> >(*this);
}




template<class R> inline
std::pair< Rational, Box<R> >
StandardIntegrator< Zonotope<R> >::
flow_bounds(const VectorField<R>& vf, 
            const Zonotope<R>& es,
            const Rational& t) const
{
  return this->_bounder->flow_bounds(vf,es.bounding_box(),t);
}


template<class R>
Zonotope<R>
StandardIntegrator< Zonotope<R> >::
integration_step(const VectorField<R>& vector_field, 
                 const Zonotope<R>& initial_set, 
                 const Rational& step_size, 
                 const Box<R>& flow_bounding_box) const
{
  
  
  
  
  uint verbosity=0;

  ARIADNE_LOG(5,"StandardIntegrator::integration_step(VectorField,Zonotope,Time,Box)\n");
  ARIADNE_LOG(6,"flow_bounding_box="<<flow_bounding_box<<"\n");
  ARIADNE_LOG(6,"initial_set="<<initial_set<<"\n");
  AffineModel<R> affine_flow_model=this->_flower->affine_flow_model(vector_field,initial_set.centre(),initial_set.bounding_box(),step_size,flow_bounding_box);
  TaylorModel<R> taylor_flow_model=this->_flower->taylor_flow_model(vector_field,initial_set.centre(),initial_set.bounding_box(),step_size,flow_bounding_box);
  ARIADNE_LOG(6,"affine_flow_model="<<affine_flow_model<<"\n");
  ARIADNE_LOG(6,"taylor_flow_model="<<taylor_flow_model<<"\n");
  ARIADNE_LOG(6,"affine_flow_model="<<affine_flow_model<<"\n");
  Zonotope<R> flow_set=apply(affine_flow_model,initial_set);
  ARIADNE_LOG(6,"flow_set="<<flow_set<<"\n");

  return flow_set;
}




template<class R>
Zonotope<R> 
StandardIntegrator< Zonotope<R> >::
reachability_step(const VectorField<R>& vector_field, 
                  const Zonotope<R>& initial_set,
                  const Rational& step_size,
                  const Box<R>& bounding_box) const
{
  
  
  
  

  uint verbosity=0;

  ARIADNE_LOG(6,"StandardIntegrator::reachability_step(VectorField,Zonotope<Interval>,Interval,Box) const\n");
  Rational half_step_size=step_size/2;

  AffineModel<R> flow_model=this->_flower->affine_flow_model(vector_field,initial_set.centre(),initial_set.bounding_box(),half_step_size,bounding_box);
  Point<I> phic=flow_model.value();
  Matrix<I> Dphi=flow_model.jacobian();
  Matrix<I> gen=Dphi*initial_set.generators();
  Vector<I> hhf=I(half_step_size)*vector_field(bounding_box);
  Vector<I> err=Dphi*(I(-1,1)*initial_set.error());

  Zonotope<R> result(phic+err,concatenate_columns(gen,hhf));

  return result;
}



template<class R>
std::ostream&
StandardIntegrator< Zonotope<R> >::write(std::ostream& os) const
{
  return os << "StandardIntegrator<Zonotope>( )";
}





}
