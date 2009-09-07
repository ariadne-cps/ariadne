/***************************************************************************
 *            solver_submodule.cc
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

#include "function_interface.h"
#include "solver_interface.h"
#include "solver.h"
#include "taylor_function.h"

#include "integrator_interface.h"
#include "integrator.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

#include "utilities.h"

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
    Vector<Interval> zero(const VectorFunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("zero")(); }
    Vector<Interval> fixed_point(const VectorFunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("fixed_point")(); }
    Set< Vector<Interval> > solve(const VectorFunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("solve")(); }
    List<VectorTaylorFunction> implicit(const VectorFunctionInterface& f, const Vector<Interval>& pd, const Vector<Interval>& bx) const {
        return this->get_override("implicit")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

class IntegratorWrapper
  : public IntegratorInterface, public wrapper< IntegratorInterface >
{
  public:
//    IntegratorInterface* clone() const {
//        return this->get_override("clone")(); }
    Pair<Float,IVector> flow_bounds(const VectorFunctionInterface&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("flow_bounds")(); }
    VectorTaylorFunction flow(const VectorFunctionInterface& vector_field,const IVector&,const Float&) const {
        return this->get_override("flow")(); }
    VectorTaylorFunction flow(const VectorFunctionInterface&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("flow")(); }
    VectorTaylorFunction time_step(const VectorFunctionInterface&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("time_step")(); }
//    std::ostream& write(std::ostream&) const {
//        return this->get_override("write")(); }
};



}


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("zero",(Vector<Interval>(SolverInterface::*)(const VectorFunctionInterface&,const Vector<Interval>&)const)&SolverInterface::zero);
    solver_wrapper_class.def("fixed_point",(Vector<Interval>(SolverInterface::*)(const VectorFunctionInterface&,const Vector<Interval>&)const)&SolverInterface::fixed_point);
    solver_wrapper_class.def("solve",(Set< Vector<Interval> >(SolverInterface::*)(const VectorFunctionInterface&,const Vector<Interval>&)const)&SolverInterface::solve);

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());

    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
    krawczyk_solver_class.def("implicit",(List<VectorTaylorFunction>(KrawczykSolver::*)(const VectorFunctionInterface&,const Vector<Interval>&,const Vector<Interval>&)const) &KrawczykSolver::implicit);

}



void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    class_<TaylorIntegrator, bases<IntegratorInterface> > taylor_integrator_class("TaylorIntegrator",init<unsigned int>());
    taylor_integrator_class.def("flow",(VectorTaylorFunction(TaylorIntegrator::*)(const VectorFunctionInterface&,const Vector<Interval>&,const Vector<Interval>&,const Float&)const)&TaylorIntegrator::flow);
}


void solver_submodule()
{
    to_python_list< Set< Vector<Interval> > >();
    to_python< List< VectorTaylorFunction > >();

    export_solver();
    export_integrator();

}


