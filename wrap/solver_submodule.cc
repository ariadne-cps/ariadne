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

#include "function.h"
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
    Vector<Interval> zero(const VectorFunction& f, const Vector<Interval>& bx) const {
        return this->get_override("zero")(); }
    Vector<Interval> fixed_point(const VectorFunction& f, const Vector<Interval>& bx) const {
        return this->get_override("fixed_point")(); }
    Vector<Interval> solve(const VectorFunction& f, const Vector<Interval>& bx) const {
        return this->get_override("solve")(); }
    VectorTaylorFunction implicit(const VectorFunction& f, const Vector<Interval>& pd, const Vector<Interval>& bx) const {
        return this->get_override("implicit")(); }
    ScalarTaylorFunction implicit(const ScalarFunction& f, const Vector<Interval>& pd, const Interval& ivl) const {
        return this->get_override("implicit")(); }
    Set< Vector<Interval> > solve_all(const VectorFunction& f, const Vector<Interval>& bx) const {
        return this->get_override("solve_all")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};


class IntegratorWrapper
  : public IntegratorInterface, public wrapper< IntegratorInterface >
{
  public:
    IntegratorInterface* clone() const {
        return this->get_override("clone")(); }
    Pair<Float,IVector> flow_bounds(const VectorFunction&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("flow_bounds")(); }
    VectorTaylorFunction flow(const VectorFunction& vector_field,const IVector&,const Float&) const {
        return this->get_override("flow")(); }
    VectorTaylorFunction flow(const VectorFunction&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("flow")(); }
    VectorTaylorFunction time_step(const VectorFunction&,const IVector&,const IVector&,const Float&) const {
        return this->get_override("time_step")(); }
//    std::ostream& write(std::ostream&) const {
//        return this->get_override("write")(); }
};



}


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",pure_virtual((Vector<Interval>(SolverInterface::*)(const VectorFunction&,const Vector<Interval>&)const) &SolverInterface::solve));
    solver_wrapper_class.def("implicit",pure_virtual((VectorTaylorFunction(SolverInterface::*)(const VectorFunction&,const Vector<Interval>&,const Vector<Interval>&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("implicit",pure_virtual((ScalarTaylorFunction(SolverInterface::*)(const ScalarFunction&,const Vector<Interval>&,const Interval&)const) &SolverInterface::implicit));
    solver_wrapper_class.def("solve_all",pure_virtual((Set< Vector<Interval> >(SolverInterface::*)(const VectorFunction&,const Vector<Interval>&)const) &SolverInterface::solve_all));
    //solver_wrapper_class.def(self_ns::str(self));

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}



void export_integrator()
{
    class_<IntegratorWrapper, boost::noncopyable> integrator_wrapper_class("IntegratorInterface");
    class_<TaylorIntegrator, bases<IntegratorInterface> > taylor_integrator_class("TaylorIntegrator",init<unsigned int>());
    taylor_integrator_class.def("flow",(VectorTaylorFunction(TaylorIntegrator::*)(const VectorFunction&,const Vector<Interval>&,const Float&)const)&TaylorIntegrator::flow);
    taylor_integrator_class.def("flow",(VectorTaylorFunction(TaylorIntegrator::*)(const VectorFunction&,const Vector<Interval>&,const Vector<Interval>&,const Float&)const)&TaylorIntegrator::flow);
}


void solver_submodule()
{
    to_python_list< Set< Vector<Interval> > >();
    to_python< List< Vector<Interval> > >();

    export_solver();
    export_integrator();

}


