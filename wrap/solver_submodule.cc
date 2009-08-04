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

#include "expression_interface.h"
#include "function_interface.h"
#include "solver_interface.h"
#include "solver.h"

#include <boost/python.hpp>
using namespace boost::python;

using namespace Ariadne;

#include "utilities.h"

namespace Ariadne {

template<class C>
struct container_to_python_list
{
    static PyObject* convert(const C& c) {
        boost::python::list result;
        for(typename C::const_iterator p=c.begin(); p!=c.end(); ++p) {
            result.append(boost::python::object(*p));
        }
        return boost::python::incref(boost::python::tuple(result).ptr());
    }
    static const PyTypeObject* get_pytype() { return &PyList_Type; }
};


class SolverWrapper
  : public SolverInterface, public wrapper< SolverInterface >
{
  public:
    SolverInterface* clone() const { return this->get_override("clone")(); }
    void set_maximum_error(double) { this->get_override("set_maximum_error")(); }
    double maximum_error() const { return this->get_override("maximum_error")(); }
    void set_maximum_number_of_steps(uint) { this->get_override("set_maximum_number_of_steps")(); }
    uint maximum_number_of_steps() const { return this->get_override("maximum_number_of_steps")(); }
    Vector<Interval> solve(const FunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("solve")(); }
    Set< Vector<Interval> > solve_all(const FunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("solve_all")(); }
    Vector<Interval> fixed_point(const FunctionInterface& f, const Vector<Interval>& bx) const {
        return this->get_override("fixed_point")(); }
    std::ostream& write(std::ostream&) const { return this->get_override("write")(); }
};

}


void export_solver()
{
    class_<SolverWrapper, boost::noncopyable> solver_wrapper_class("SolverInterface");
    solver_wrapper_class.def("solve",&SolverInterface::solve);
    solver_wrapper_class.def("solve_all",&SolverInterface::solve_all);

    class_<IntervalNewtonSolver, bases<SolverInterface> > interval_newton_solver_class("IntervalNewtonSolver",init<double,unsigned int>());
    class_<KrawczykSolver, bases<SolverInterface> > krawczyk_solver_class("KrawczykSolver",init<double,unsigned int>());
}




void solver_submodule()
{
    boost::python::to_python_converter< Set< Vector<Interval> >, set_to_python_list< Vector<Interval> > >();
    export_solver();
}


