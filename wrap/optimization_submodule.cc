/***************************************************************************
 *            optimization_submodule.cc
 *
 *  Copyright 2008  Pieter Collins
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
 
#include "config.h"

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

#include "linear_programming.h"

#include "utilities.h"

#include <boost/python.hpp>

using namespace boost::python;

using namespace Ariadne;


template<class X>
boost::python::tuple
python_compute_basis(const Matrix<X>& A) {
    array<size_t> p;
    Matrix<X> B;
    make_lpair(p,B)=compute_basis(A);
    boost::python::list l;
    for(size_t i=0; i!=p.size(); ++i) {
        l.append(p[i]);
    }
    return boost::python::make_tuple(l,B);
}

template<class T> T get(const array<T>& ary, size_t i) { return ary[i]; }
template<class T> void set(array<T>& ary, size_t i, const T& t) { ary[i]=t; }

template<class T>
void export_internal_array(const char* name)
{
    class_< array<T> > array_class(name,no_init);
    array_class.def("__len__", &array<T>::size);
    array_class.def("__getitem__",&get<T>);
    array_class.def(boost::python::self_ns::str(self));
}


void export_variable_type()
{
    typedef array<VariableType> VariableTypeArray;
    
    enum_<VariableType> variable_enum("VariableType");
    variable_enum.value("BASIS", BASIS);
    variable_enum.value("LOWER", LOWER);
    variable_enum.value("UPPER", UPPER);
}

template<class X>
void export_linear_programming()
{
    typedef array<size_t> SizeArray;

    to_python< std::pair< array<size_t>, Matrix<X> > >();

    def("lpstep",(bool(*)(const Matrix<X>&,const Vector<X>&,const Vector<X>&,SizeArray&,Matrix<X>&,Vector<X>&)) &lpstep);
    def("lpstep",(bool(*)(const Matrix<X>&,const Vector<X>&,const Vector<X>&,const Vector<X>&,const Vector<X>&,array<VariableType>&,array<size_t>&,Matrix<X>&,Vector<X>&)) &lpstep);
  

    def("primal_feasible",(tribool(*)(const Matrix<X>&,const Vector<X>&)) &primal_feasible);
    def("dual_feasible",(tribool(*)(const Matrix<X>&,const Vector<X>&)) &dual_feasible);
    def("constrained_feasible",(tribool(*)(const Matrix<X>&,const Vector<X>&,const Vector<X>&,const Vector<X>&)) &constrained_feasible);

    def("verify_primal_feasibility",(tribool(*)(const Matrix<X>&,const Vector<X>&,const array<VariableType>&)) &verify_primal_feasibility);
    def("verify_dual_feasibility",(tribool(*)(const Matrix<X>&,const Vector<X>&,const array<VariableType>&)) &verify_dual_feasibility);
    def("verify_constrained_feasibility",(tribool(*)(const Matrix<X>&,const Vector<X>&,const Vector<X>&,const Vector<X>&,const array<VariableType>&)) &verify_constrained_feasibility);

    def("compute_basis",(std::pair< SizeArray, Matrix<X> >(*)(const Matrix<X>&)) &compute_basis<X>);

    def("constrained_feasible_by_enumeration",(tribool(*)(const Matrix<X>&,const Vector<X>&,const Vector<X>&,const Vector<X>&)) &constrained_feasible_by_enumeration);
}
   
#ifdef HAVE_GMPXX_H
template void export_linear_programming<Rational>();
#endif


void optimization_submodule() {
    export_variable_type();
    export_array<size_t>("SizeArray");
    export_internal_array<VariableType>("VariableTypeArray");
    export_linear_programming<Float>();
#ifdef HAVE_GMPXX_H
    export_linear_programming<Rational>();
#endif
}
