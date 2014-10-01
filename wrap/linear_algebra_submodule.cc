/***************************************************************************
 *            linear_algebra_submodule.cc
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

#include "boost_python.h"
#include "utilities.h"

#include "config.h"

#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include "numeric.h"
#include "vector.h"
#include "matrix.h"

using namespace boost::python;
using namespace Ariadne;


namespace Ariadne {


template<class X>
X __vgetitem__(const Vector<X>& v, int i)
{
    if(i<0) { i+=v.size(); }
    ARIADNE_ASSERT_MSG(0<=i && uint(i)<v.size(),"v="<<v<<" i="<<i);
    return v[i];
}


template<class X>
Vector<X> __vgetslice__(const Vector<X>& v, int start, int stop)
{
    if(start<0) { start+=v.size(); }
    if(stop<0) { stop+=v.size(); }
    ARIADNE_ASSERT(0<=start && start<=stop && uint(stop)<=v.size());
    return project(v,range(start,stop));
}


template<class X>
void __vsetitem__(Vector<X>& v, int i, const X& x)
{
    if(i<0) { i+=v.size(); }
    ARIADNE_ASSERT(0<=i && uint(i)<v.size());
    v[i]=x;
}


template<class X>
X __mgetitem__(const Matrix<X>& A, const boost::python::tuple& tup)
{
    uint i=boost::python::extract<uint>(tup[0]);
    uint j=boost::python::extract<uint>(tup[1]);
    return A[i][j];
}

template<class X>
void __msetitem__(Matrix<X>& A, const boost::python::tuple& tup, const X& x)
{
    uint i=boost::python::extract<uint>(tup[0]);
    uint j=boost::python::extract<uint>(tup[1]);
    A[i][j]=x;
}

template<class X> Matrix<X> transpose(const Matrix<X>& A) { return MatrixTranspose< Matrix<X> >(A); }


template<class X>
struct from_python< Vector<X> >
{
    from_python() {
        converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >());
    }

    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data)
    {
        list lst=extract<list>(obj_ptr);
        void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*)   data)->storage.bytes;
        Vector<X> res(len(lst));
        for(uint i=0; i!=res.size(); ++i) { res[i]=extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
        data->convertible = storage;
    }
};


template<class X>
struct from_python< Matrix<X> >
{
    from_python() {
        converter::registry::push_back(&convertible,&construct,type_id< Matrix<X> >());
    }

    static void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data)
    {
        list rows=extract<list>(obj_ptr);
        Matrix<X> res(len(rows),len(extract<list>(rows[0])));
        for(uint i=0; i!=res.row_size(); ++i) {
            list elmnts=extract<list>(rows[i]);
            ARIADNE_ASSERT(uint(len(elmnts))==res.column_size());
            for(uint j=0; j!=res.column_size(); ++j) {
                res[i][j]=extract<X>(elmnts[j]); } }
        void* storage = ((converter::rvalue_from_python_storage< Matrix<X> >*)   data)->storage.bytes;
        new (storage) Matrix<X>(res);
        data->convertible = storage;
    }
};

std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ApproximateFloat>& repr);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ValidatedFloat>& repr);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<ExactFloat>& repr);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<RawFloat>& repr);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Rational>& repr);
std::ostream& operator<<(std::ostream& os, const PythonRepresentation<Real>& repr);

template<class X> std::ostream& operator<<(std::ostream& os, const PythonRepresentation< Vector<X> >& repr) {
    Vector<X> const& v=repr.reference();
    os << "["; for(uint i=0; i!=v.size(); ++i) { if(i!=0) { os << ","; } os << python_representation(v[i]); } os << "]"; return os;
}

template<class X> std::ostream& operator<<(std::ostream& os, const PythonRepresentation< Matrix<X> >& repr) {
    Matrix<X> const& A=repr.reference();
    os << "["; for(uint i=0; i!=A.row_size(); ++i) { if(i!=0) { os << ","; } os << "["; for(uint j=0; j!=A.column_size(); ++j) {
        if(j!=0) { os << ","; } os << python_representation(A[i][j]); } os << "]"; } os << "]"; return os;
}

template<class X> auto __norm__(const Vector<X>& v) -> decltype(norm(v)) { return norm(v); }
template<class X> auto __dot__(const Vector<X>& v1, const Vector<X>& v2) -> decltype(dot(v1,v2)) { return dot(v1,v2); }
template<class X> Vector<X> __join__(const Vector<X>& v1, const Vector<X>& v2) { return join(v1,v2); }
template<class X> Vector<X> __join__(const Vector<X>& v1, const X& s2) { return join(v1,s2); }
template<class X> Vector<X> __join__(const X& s1, const Vector<X>& v2) { return join(s1,v2); }
template<class X> Vector<X> __join__(const X& s1, const X& s2) { return join(s1,s2); }

} // namespace Ariadne


template<class X>
void export_vector_class(class_<Vector<X> >& vector_class)
{
    vector_class.def(init<int>());
    vector_class.def(init<int,X>());
    vector_class.def("size", &Vector<X>::size);
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", &__vsetitem__<X>);
    vector_class.def("__setitem__", &__vsetitem__<X>);
    vector_class.def("__getitem__", &__vgetitem__<X>);
    vector_class.def("__getslice__", &__vgetslice__<X>);
    vector_class.def("__eq__", &__eq__<bool,Vector<X>,Vector<X> >);
    vector_class.def("__ne__", &__ne__<bool,Vector<X>,Vector<X> >);
    vector_class.def("__pos__", &__pos__< Vector<X>, Vector<X> >);
    vector_class.def("__neg__", &__neg__< Vector<X>, Vector<X> >);
    vector_class.def("__str__",&__cstr__< Vector<X> >);
    vector_class.def("__repr__",&__repr__< Vector<X> >);
    vector_class.def("unit",&Vector<X>::unit);
    vector_class.def("basis",&Vector<X>::basis);
    vector_class.staticmethod("unit");
    vector_class.staticmethod("basis");

    def("join",(Vector<X>(*)(const Vector<X>&,const Vector<X>&)) &__join__);
    def("join",(Vector<X>(*)(const Vector<X>&,const X&)) &__join__);
    def("join",(Vector<X>(*)(const X&,const Vector<X>&)) &__join__);
    def("join",(Vector<X>(*)(const X&,const X&)) &__join__);

    from_python< Vector<X> >();
    to_python< Array< Vector<X> > >();


}

template<class X, class Y>
void export_vector_conversion(class_<Vector<X> >& vector_class)
{
    //vector_class.def(init<int,Y>());
    vector_class.def(init< Vector<Y> >());
}

template<class R, class X, class Y>
void export_vector_arithmetic(class_<Vector<X> >& vector_class)
{
    vector_class.def("__add__",__add__< Vector<R>, Vector<X>, Vector<Y> >);
    vector_class.def("__sub__",__sub__< Vector<R>, Vector<X>, Vector<Y> >);
    vector_class.def("__rmul__",__rmul__< Vector<R>, Vector<X>, Y >);
    vector_class.def("__mul__",__mul__< Vector<R>, Vector<X>, Y >);
    vector_class.def("__div__",__div__< Vector<R>, Vector<X>, Y >);
    // Don't use boost::python style operators self+other (as below) because of
    // below) because of need to convert result expressions to Vector<R>.
    // vector_class.def(self+=Vector<Y>());
    // vector_class.def(self*=Y());
}


template<class X> void export_vector()
{
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<X> >());
    export_vector_class<X>(vector_class);
    export_vector_conversion<X,X>(vector_class);
    export_vector_arithmetic<X,X,X>(vector_class);
}


template<> void export_vector<ExactFloat>()
{
    class_< Vector<ExactFloat> > exact_vector_class("ExactFloatVector",init< Vector<ExactFloat> >());
    export_vector_class<ExactFloat>(exact_vector_class);
}

template<> void export_vector<ValidatedFloat>()
{
    class_< Vector<ValidatedFloat> > validated_vector_class("ValidatedFloatVector",init< Vector<ValidatedFloat> >());
    export_vector_class<ValidatedFloat>(validated_vector_class);
    export_vector_conversion<ValidatedFloat,ExactFloat>(validated_vector_class);
    export_vector_arithmetic<ValidatedFloat,ValidatedFloat,ValidatedFloat>(validated_vector_class);
    def("norm",&__norm__<ValidatedFloat>);
    def("dot", &__dot__<ValidatedFloat>);

    implicitly_convertible< Vector<ExactFloat>, Vector<ValidatedFloat> >();
}

template<> void export_vector<ApproximateFloat>()
{
    class_< Vector<ApproximateFloat> > approximate_vector_class("FloatVector",init< Vector<ApproximateFloat> >());
    export_vector_class<ApproximateFloat>(approximate_vector_class);
    export_vector_conversion<ApproximateFloat,ValidatedFloat>(approximate_vector_class);
    export_vector_arithmetic<ApproximateFloat,ApproximateFloat,ApproximateFloat>(approximate_vector_class);
    //export_vector_arithmetic<ValidatedFloat,ApproximateFloat,ValidatedFloat>(approximate_vector_class);
    approximate_vector_class.def("__rmul__",__rmul__< Vector<ApproximateFloat>, Vector<ApproximateFloat>, ApproximateFloat >);
    approximate_vector_class.def("__mul__",__mul__< Vector<ApproximateFloat>, Vector<ApproximateFloat>, ApproximateFloat >);
    approximate_vector_class.def("__div__",__div__< Vector<ApproximateFloat>, Vector<ApproximateFloat>, ApproximateFloat >);
}




template<class X>
void export_matrix_class(class_<Matrix<X> >& matrix_class)
{
    typedef uint SizeType;

    matrix_class.def(init<int,int>());
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", &__msetitem__<X>);
    matrix_class.def("__getitem__", &__mgetitem__<X>);
    matrix_class.def("__pos__", &__pos__< Matrix<X>, Matrix<X> >);
    matrix_class.def("__neg__", &__neg__< Matrix<X>, Matrix<X> >);
    matrix_class.def("__str__",&__cstr__< Matrix<X> >);
    matrix_class.def("__repr__",&__repr__<Matrix<X> >);

    matrix_class.def("identity",(Matrix<X>(*)(size_t)) &Matrix<X>::identity);
    matrix_class.staticmethod("identity");

    from_python< Matrix<X> >();
}


template<class X, class Y>
void export_matrix_conversion(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(init< Matrix<Y> >());
}

template<class R, class X, class Y>
void export_matrix_arithmetic(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def("__add__", &__add__< Matrix<R>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__sub__", &__sub__< Matrix<R>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__mul__", &__mul__< Matrix<R>, Matrix<X>, Y >);
    matrix_class.def("__rmul__", &__rmul__< Matrix<R>, Matrix<X>, Y >);
    matrix_class.def("__div__", &__div__< Matrix<R>, Matrix<X>, Y >);
    matrix_class.def("__mul__", &__mul__< Vector<R>, Matrix<X>, Vector<Y> >);
    matrix_class.def("__mul__", &__mul__< Matrix<R>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__rmul__", &__rmul__< Vector<R>, Matrix<X>, Vector<Y> >);
}

template<class X>
void export_matrix_operations(class_<Matrix<X> >& matrix_class)
{
    def("norm",(X(*)(const Matrix<X>&)) &norm);
    def("transpose",(Matrix<X>(*)(const Matrix<X>&)) &transpose);

    def("inverse",(Matrix<X>(*)(const Matrix<X>&)) &inverse);
    def("solve",(Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &solve);
    def("solve",(Vector<X>(*)(const Matrix<X>&,const Vector<X>&)) &solve);
}


template<class X> void export_matrix()
{
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init< Matrix<X> >());
    export_matrix_class<X>(matrix_class);
    export_matrix_conversion<X,X>(matrix_class);
    export_matrix_arithmetic<X,X,X>(matrix_class);
    export_matrix_arithmetic<X,X,double>(matrix_class);
    export_matrix_operations<X>(matrix_class);


}

template<> void export_matrix<ExactFloat>()
{
    class_< Matrix<ExactFloat> > matrix_class(python_name<ExactFloat>("Matrix"),no_init);
    matrix_class.def(init<PivotMatrix>());
    export_matrix_class<ExactFloat>(matrix_class);
}

template<> void export_matrix<ValidatedFloat>()
{
    class_< Matrix<ValidatedFloat> > matrix_class(python_name<ValidatedFloat>("Matrix"),no_init);
    export_matrix_class<ValidatedFloat>(matrix_class);
    export_matrix_conversion<ValidatedFloat,ExactFloat>(matrix_class);
    export_matrix_arithmetic<ValidatedFloat,ValidatedFloat,ValidatedFloat>(matrix_class);
    //export_matrix_arithmetic<ValidatedFloat,ValidatedFloat,ApproximateFloat>(matrix_class);
    def("gs_inverse", (Matrix<ValidatedFloat>(*)(const Matrix<ValidatedFloat>&)) &gs_inverse);
    def("lu_inverse", (Matrix<ValidatedFloat>(*)(const Matrix<ValidatedFloat>&)) &lu_inverse);
    def("gs_solve", (Vector<ValidatedFloat>(*)(const Matrix<ValidatedFloat>&,const Vector<ValidatedFloat>&)) &gs_solve);
    def("gs_solve", (Matrix<ValidatedFloat>(*)(const Matrix<ValidatedFloat>&,const Matrix<ValidatedFloat>&)) &gs_solve);
    def("lu_solve", (Matrix<ValidatedFloat>(*)(const Matrix<ValidatedFloat>&,const Matrix<ValidatedFloat>&)) &lu_solve);

    //implicitly_convertible< Matrix<ApproximateFloat>, Matrix<ValidatedFloat> >();
}

template<> void export_matrix<ApproximateFloat>()
{
    class_< Matrix<ApproximateFloat> > matrix_class(python_name<ApproximateFloat>("Matrix"),no_init);
    export_matrix_class<ApproximateFloat>(matrix_class);
    export_matrix_conversion<ApproximateFloat,ValidatedFloat>(matrix_class);
    export_matrix_arithmetic<ApproximateFloat,ApproximateFloat,ApproximateFloat>(matrix_class);
    //export_matrix_arithmetic<ValidatedFloat,ApproximateFloat,ValidatedFloat>(matrix_class);

    def("triangular_decomposition",&triangular_decomposition);
    def("orthogonal_decomposition", &orthogonal_decomposition);
    def("triangular_factor",&triangular_factor);
    def("triangular_multiplier", &triangular_multiplier);
    def("row_norms",(Vector<ApproximateFloat>(*)(const Matrix<ApproximateFloat>&)) &row_norms);
    def("normalise_rows",(Matrix<ApproximateFloat>(*)(const Matrix<ApproximateFloat>&)) &normalise_rows);

    to_python< Tuple<FloatMatrix,FloatMatrix,PivotMatrix> >();
}


template<class X> void export_diagonal_matrix()
{
    class_< DiagonalMatrix<X> > diagonal_matrix_class(python_name<X>("DiagonalMatrix"),no_init);
    diagonal_matrix_class.def(init< uint >());
    diagonal_matrix_class.def(init< Vector<X> >());
    diagonal_matrix_class.def("__setitem__", &DiagonalMatrix<X>::set);
    diagonal_matrix_class.def("__getitem__", &DiagonalMatrix<X>::get, return_value_policy<copy_const_reference>());
    diagonal_matrix_class.def("__add__", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&,const DiagonalMatrix<X>&)) operator+ );
    diagonal_matrix_class.def("__sub__", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&,const DiagonalMatrix<X>&)) operator- );
    diagonal_matrix_class.def("__mul__", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&,const DiagonalMatrix<X>&)) operator* );
    diagonal_matrix_class.def("__div__", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&,const DiagonalMatrix<X>&)) operator/ );
    diagonal_matrix_class.def("__mul__", (Vector<X>(*)(const DiagonalMatrix<X>&,const Vector<X>&)) operator* );
    diagonal_matrix_class.def("__mul__", (Matrix<X>(*)(const DiagonalMatrix<X>&,const Matrix<X>&)) operator* );
    diagonal_matrix_class.def("__str__",&__cstr__< DiagonalMatrix<X> >);
    //diagonal_matrix_class.def("__mul__", &__mul__< DiagonalMatrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    //diagonal_matrix_class.def("__div__", &__div__< DiagonalMatrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    //diagonal_matrix_class.def("__mul__", &__mul__< Vector<X>, DiagonalMatrix<X>, Vector<X> >);
    //diagonal_matrix_class.def("__mul__", &__mul__< Matrix<X>, DiagonalMatrix<X>, Matrix<X> >);
    //diagonal_matrix_class.def("__rmul__", &__prod__< Vector<X>, Vector<X>, DiagonalMatrix<X> >);
    //diagonal_matrix_class.def("__rmul__", &__prod__< Matrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    //def("inverse", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&)) &inverse<X>);
}

void export_pivot_matrix()
{

    implicitly_convertible< PivotMatrix, Matrix<ExactFloat> >();

    class_<PivotMatrix> pivot_matrix_class("PivotMatrix",no_init);
    pivot_matrix_class.def("__str__",&__cstr__<PivotMatrix>);
    pivot_matrix_class.def("__repr__",&__cstr__<PivotMatrix>);
}


template void export_vector<ApproximateFloat>();
template void export_vector<ValidatedFloat>();
template void export_matrix<ApproximateFloat>();
template void export_matrix<ValidatedFloat>();

template void export_diagonal_matrix<ApproximateFloat>();

#ifdef HAVE_GMPXX_H
template void export_vector<Rational>();
template void export_matrix<Rational>();
#endif


void linear_algebra_submodule() {
    export_vector<ApproximateFloat>();
    export_vector<ValidatedFloat>();

    export_matrix<ApproximateFloat>();
    export_matrix<ValidatedFloat>();

    export_pivot_matrix();
    export_diagonal_matrix<ApproximateFloat>();

#ifdef HAVE_GMPXX_H
    export_vector<Rational>();
    export_matrix<Rational>();
#endif
}
