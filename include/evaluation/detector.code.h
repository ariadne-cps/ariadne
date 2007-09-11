/***************************************************************************
 *            detector.code.h
 *
 *  Copyright  2006  Alberto Casagrande, Pieter Collins
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
 
#include "system/vector_field_interface.h"
#include "geometry/constraint_interface.h"
#include "geometry/basic_set_interface.h"

#include "geometry/point.h"
#include "geometry/rectangle.h"
#include "geometry/zonotope.h"
#include "geometry/constraint.h"

#include "evaluation/exceptions.h"
#include "../evaluation/time_model.h"
#include "evaluation/bounder.h"
#include "evaluation/lohner_integrator.h"

#include "output/logging.h"

namespace Ariadne {

namespace Evaluation { static int& verbosity = detector_verbosity; }

template<class R>
Evaluation::Detector<R>::~Detector()
{
}

template<class R>
Evaluation::Detector<R>::Detector()
{
}

template<class R>
Evaluation::Detector<R>::Detector(const Detector<R>& det) 
{
}

template<class R>
Evaluation::Detector<R>*
Evaluation::Detector<R>::clone() const
{
  return new Detector<R>(*this);
}




template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::value(const Geometry::ConstraintInterface<R>& c, 
                               const Geometry::Rectangle<R>& r) const
{
  Geometry::Point<I> pt=r;
  return c.value(pt);
}


template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::value(const Geometry::ConstraintInterface<R>& c, 
                               const Geometry::Zonotope<R,R>& z) const
{
  const Geometry::DifferentiableConstraintInterface<R>& dc=dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c);
  if(&dc) {
    LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
    const Geometry::Point<R>& zc=z.centre();
    const LinearAlgebra::Matrix<R>& zG=z.generators();
    Geometry::Rectangle<R> bb=z.bounding_box();
    Geometry::Point<I> bpt(bb);
    Numeric::Interval<R> v=dc.value(zc)+LinearAlgebra::inner_product(dc.gradient(bb)*zG,e);
    return v;
  } else {
    return this->value(c,z.bounding_box());
  }
}

template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::value(const Geometry::ConstraintInterface<R>& c, 
                               const Geometry::Zonotope<I,R>& z) const
{
  const Geometry::DifferentiableConstraintInterface<R>& dc=dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c);
  if(&dc) {
    LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
    const Geometry::Point<I>& zc=z.centre();
    const LinearAlgebra::Matrix<R>& zG=z.generators();
    Geometry::Rectangle<R> bb=z.bounding_box();
    Geometry::Point<I> bpt(bb);
    Numeric::Interval<R> v=dc.value(zc)+LinearAlgebra::inner_product(dc.gradient(bb)*zG,e);
    return v;
  } else {
    return this->value(c,z.bounding_box());
  }
}

template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::value(const Geometry::ConstraintInterface<R>& c, 
                               const Geometry::Zonotope<I,I>& z) const
{
  const Geometry::DifferentiableConstraintInterface<R>& dc=dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c);
  if(&dc) {
    LinearAlgebra::Vector<I> e(z.number_of_generators(),I(-1,1));
    const Geometry::Point<I>& zc=z.centre();
    const LinearAlgebra::Matrix<I>& zG=z.generators();
    Geometry::Rectangle<R> bb=z.bounding_box();
    Geometry::Point<I> bpt(bb);
    Numeric::Interval<R> v=dc.value(zc)+LinearAlgebra::inner_product(dc.gradient(bb)*zG,e);
    return v;
  } else {
    return this->value(c,z.bounding_box());
  }
}



template<class R>
tribool 
Evaluation::Detector<R>::forces(const Geometry::ConstraintInterface<R>& c1,
                                const Geometry::ConstraintInterface<R>& c2,
                                const Geometry::Rectangle<R>& dom) const
{
  const Geometry::DifferentiableConstraintInterface<R>& dc1=dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c1);
  const Geometry::DifferentiableConstraintInterface<R>& dc2=dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c2);
  
  if(&dc1 && &dc2) {
    Geometry::Point<I> centre = dom.centre();
    Geometry::Point<I> bounding_point = dom;
    LinearAlgebra::Vector<I> vector_radius = bounding_point-centre;
    Numeric::Interval<R> centre_difference = dc1.value(centre) - dc2.value(centre);
    LinearAlgebra::Vector<I> gradient_difference = dc1.gradient(bounding_point)-dc2.gradient(bounding_point);
    Numeric::Interval<R> interval_difference = centre_difference+inner_product(gradient_difference,vector_radius);
    if(interval_difference<0) { return true; }
    else if ( interval_difference>0) { return false; }
    else { return indeterminate; }
  } else {
    return this->value(c1,dom) < this->value(c2,dom);
  }
}


template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::normal_derivative(const System::VectorFieldInterface<R>& vf, 
                                           const Geometry::DifferentiableConstraintInterface<R>& dc, 
                                           const Geometry::Point<I>& pt) const
{  
  return LinearAlgebra::inner_product(dc.gradient(pt),vf(pt));
}


template<class R>
Numeric::Interval<R> 
Evaluation::Detector<R>::crossing_time(const System::VectorFieldInterface<R>& vf, 
                                       const Geometry::ConstraintInterface<R>& c, 
                                       const Geometry::Point<I>& pt,
                                       const Geometry::Rectangle<R>& bb) const
{  
  const Geometry::DifferentiableConstraintInterface<R>& dc=
    dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c);
  assert(&dc);
  Geometry::Point<I> bpt=bb;
  I icv = c.value(pt);
  I nd = LinearAlgebra::inner_product(dc.gradient(bpt),vf(pt));
  return -icv/nd;
}


template<class R>
Evaluation::TimeModel<R> 
Evaluation::Detector<R>::crossing_time(const System::VectorFieldInterface<R>& vf, 
                                       const Geometry::ConstraintInterface<R>& c, 
                                       const Geometry::Rectangle<R>& d,
                                       const Geometry::Rectangle<R>& bb) const
{
  ARIADNE_LOG(8,"    Detector::crossing_time(...)\n");
  static const int number_of_newton_steps=2;

  const System::VectorFieldInterface<R>& dynamic=vf;
  const Geometry::DifferentiableConstraintInterface<R>& constraint=
    dynamic_cast<const Geometry::DifferentiableConstraintInterface<R>&>(c);
  const Geometry::Rectangle<R>& domain=d;
  const Geometry::Rectangle<R>& bounding_box=bb;
  Evaluation::LohnerIntegrator<R> lohner_integrator;
  Evaluation::Bounder<R> bounder;
  
  if(!&constraint) {
    throw std::runtime_error("Detector::crossing_time(...): Can only compute crossing time for differentiable constraint");
  }


  Numeric::Rational flow_time = 0;
  Geometry::Point<I> centre = domain.centre();
  LinearAlgebra::Vector<I> centre_flow_direction;
  LinearAlgebra::Vector<I> centre_constraint_gradient;
  Numeric::Interval<R> centre_normal_derivative;
  Numeric::Interval<R> centre_constraint_value;
  Numeric::Interval<R> centre_time_step;
  Geometry::Rectangle<R> centre_bounding_box;

  try {
    // Estimate crossing time for centre by taking Newton iterations
    for(int i=0; i!=number_of_newton_steps; ++i) {
      centre = Geometry::midpoint(centre);
      centre_flow_direction = dynamic(centre);
      centre_constraint_value = constraint.value(centre);
      centre_constraint_gradient = constraint.gradient(centre);
      centre_normal_derivative = inner_product(centre_flow_direction,centre_constraint_gradient);
      if(possibly(centre_normal_derivative==0)) {
        throw NonTransverseCrossingException();
      }
      flow_time += Numeric::Rational(centre_time_step.midpoint());
      Numeric::Interval<R> time_interval=Numeric::Rational(flow_time);

      centre_bounding_box=bounder.refine_flow_bounds(dynamic,centre,bounding_box,centre_time_step.upper());
      centre_bounding_box=bounder.refine_flow_bounds(dynamic,centre,centre_bounding_box,centre_time_step.upper());
      
      centre=lohner_integrator.flow_step(dynamic,centre,flow_time,centre_bounding_box);
    }
  
    ARIADNE_LOG(9,"    estimated_centre_crossing_time="<<flow_time<<"\n");
    
    // Perform integration to close to centre
    LinearAlgebra::Vector<I> flow_direction = vf(bounding_box);
    ARIADNE_LOG(9,"    flow_direction="<<flow_direction<<"\n");
    Numeric::Interval<R>  constraint_value = constraint.value(centre);
    ARIADNE_LOG(9,"    constraint_value="<<constraint_value<<"\n");
    LinearAlgebra::Vector<I> constraint_gradient = constraint.gradient(bounding_box);
    ARIADNE_LOG(9,"    constraint_gradient="<<constraint_gradient<<"\n");
    Numeric::Interval<R> normal_derivative = -inner_product(constraint_gradient,flow_direction);
    ARIADNE_LOG(9,"    normal_derivative="<<normal_derivative<<"\n");

    // Estimate centre crossing time
    Numeric::Interval<R> centre_normal_derivative = normal_derivative;
    Numeric::Interval<R> centre_crossing_time = constraint.value(centre)/centre_normal_derivative;
    
    // Compute the gradient of the crossing times
    LinearAlgebra::Vector<I> spacial_time_gradient = constraint_gradient/normal_derivative;
    
    // Log the crossing time step and return
    TimeModel<R> spacial_crossing_time_step(centre_crossing_time, spacial_time_gradient);
    ARIADNE_LOG(9,"    spacial_crossing_time_step="<<spacial_crossing_time_step<<"\n");
    return spacial_crossing_time_step;
  }

  // Use the following code for non-transverse crossings
  catch(NonTransverseCrossingException) {
    ARIADNE_LOG(9,"   Non-transverse crossing");
    LinearAlgebra::Vector<I> flow_direction = dynamic(bounding_box);
    ARIADNE_LOG(9,"    flow_direction="<<flow_direction<<"\n");
    Numeric::Interval<R> constraint_value = constraint.value(centre);
    ARIADNE_LOG(9,"    constraint_value="<<constraint_value<<"\n");
    LinearAlgebra::Vector<I> constraint_gradient = constraint.gradient(bounding_box);
    ARIADNE_LOG(9,"    constraint_gradient="<<constraint_gradient<<"\n");
    Numeric::Interval<R> normal_derivative = -inner_product(constraint_gradient,flow_direction);
    ARIADNE_LOG(9,"    normal_derivative="<<normal_derivative<<"\n");
    
    R minimum_crossing_time = (constraint.value(bounding_box)/normal_derivative).lower();
    R maximum_crossing_time = (R(1)/R(0)).upper(); // Should be inf (infinity)
    Numeric::Interval<R> crossing_time(minimum_crossing_time,maximum_crossing_time);
    LinearAlgebra::Vector<I> spacial_time_gradient(domain.dimension());

    TimeModel<R> spacial_crossing_time_step(crossing_time, spacial_time_gradient);
    ARIADNE_LOG(9,"    spacial_crossing_time_step="<<spacial_crossing_time_step<<"\n");
    return spacial_crossing_time_step;
  } 
}



}
