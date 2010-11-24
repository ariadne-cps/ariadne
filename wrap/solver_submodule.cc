/***************************************************************************
 *      solver_submodule.cc
 *
 *  Copyright 2009  Pieter Collins
 *
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

#include <boost/python.hpp>

#include "function.h"
#include "solver_interface.h"
#include "solver.h"
#include "taylor_function.h"

#include "integrator_interface.h"
#include "integrator.h"

#include "utilities.h"

using namespace boost::python;
using namespace Ariadne;

namespace Ariadne {


class SolverWrapper
  : public SolverInterface, public wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    void set_maximum_error(double) { this->get_override("set_maximum_error")(); }
    double maximum_error() const { return this->get_override("maximum_error")(); }
    void set_maximum_number_of_steps(uint) { this->get_override("set_maximum_number_of_steps")(); }
    uint maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    IntervalVector zero(const RealVectorFunction& f, const IntervalVector& bx) const {
     return this->get_override("zero")(); }
    IntervalVector fixed_point(const RealVectorFunction& f, const IntervalVector& bx) const {
     return this->get_override("fixed_point")(); }
    IntervalVector solve(const RealVectorFunction& f, const IntervalVector& bx) const {
     return this->get_override("solve")(); }
    VectorTaylorFunction implicit(const RealVectorFunction& f, const IntervalVector& pd, const IntervalVector& bx) const {
     return this->get_override("implicit")(); }
    ScalarTaylorFunction implicit(const RealScalarFunction& f, const IntervalVector& pd, const Interval& ivl) const {
     return this->get_override("implicit")(); }
    Set< IntervalVector > solve_all(const RealVectorFunction& f, const IntervalVector& bx) const {
     return this->get_override("solve_all")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class IntegratorWrapper
  : public IntegratorInterface, public wrapper< IntegratorInterface >
{
  public:
    IntegratorInterface* clone() const {
     return this->get_override("clone")(); }
    void set_temporal_order(uint) {
     this->get_override("set_temporal_order")(); }
    void set_maximum_error(double) {
     this->get_override("set_maximum_error")(); }
    Pair<Float,IntervalVector> flow_bounds(const RealVectorFunction&,const IntervalVector&,const Float&) const {
     return this->get_override("flow_bounds")(); }
    VectorTaylorFunction flow_step(const RealVectorFunction&,const IntervalVector&,const Float&) const {
     return this->get_override("flow_step")(); }
    VectorTaylorFunction flow_step(const RealVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&) const {
     return this->get_override("flow_step")(); }
    VectorTaylorFunction flow(const RealVectorFunction& vector_field,const IntervalVector&,const Real&) const {
     return this->get_override("flow")(); }
    VectorTaylorFunction flow(const RealVectorFunction&,const IntervalVector&,const Interval&) const {
     return this->get_override("flow")(); }
//    std::ostream& write(std::ostream&) const {
//     return this->get_override("write")(); }
};


} // namespace Ariadne


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",pure_virtual((IntervalVector(SolverInterface::*)(const RealVectorFunction&,const IntervalVector&)const) &SolverInterface::solve));
    solver_wrapper_class.def("implicit",pure_virtual((VectorTaylorFunction(SolverInterface::*)(const RealVectorFunction&,const IntervalVector&,const IntervalVector&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("implicit",pure_virtual((ScalarTaylorFunction(SolverInterface::*)(const RealScalarFunction&,const IntervalVector&,const Interval&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("solve_all",pure_virtual((Set< IntervalVector >(SolverInterface::*)(const RealVectorFunction&,const IntervalVector&)const) &SolverInterface::solve_all));
    //solver_wrapper_class.def(self_ns::str(self));

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}



void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    integrator_wrapper_class.def("flow",(VectorTaylorFunction(IntegratorInterface::*)(const RealVectorFunction&,const IntervalVector&,const Real&)const)&IntegratorInterface::flow);
    integrator_wrapper_class.def("flow",(VectorTaylorFunction(IntegratorInterface::*)(const RealVectorFunction&,const IntervalVector&,const Interval&)const)&IntegratorInterface::flow);
    integrator_wrapper_class.def("flow_bounds",(Pair<Float,IntervalVector>(IntegratorInterface::*)(const RealVectorFunction&,const IntervalVector&,const Float&)const)&IntegratorInterface::flow_bounds);
    integrator_wrapper_class.def("flow_step",(VectorTaylorFunction(IntegratorInterface::*)(const RealVectorFunction&,const IntervalVector&,const Float&)const)&IntegratorInterface::flow_step);
    integrator_wrapper_class.def("flow_step",(VectorTaylorFunction(IntegratorInterface::*)(const RealVectorFunction&,const IntervalVector&,const Float&,const IntervalVector&)const)&IntegratorInterface::flow_step);
    class_<TaylorIntegrator, bases<IntegratorInterface> > taylor_integrator_class("TaylorIntegrator",init<unsigned int,double>());
}


void solver_submodule()
{
    to_python_list< Set< IntervalVector > >();
    to_python< List< IntervalVector > >();

    export_solver();
    export_integrator();

}


