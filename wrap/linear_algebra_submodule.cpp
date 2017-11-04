/***************************************************************************
 *            linear_algebra_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
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

#include "boost_python.hpp"
#include "utilities.hpp"

#include "config.h"

#include <boost/python.hpp>
#include <boost/python/slice.hpp>

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"

using namespace boost::python;
using namespace Ariadne;


namespace Ariadne {


template<class X>
X __vgetitem__(const Vector<X>& v, Int i)
{
    if(i<0) { i+=v.size(); }
    ARIADNE_ASSERT_MSG(0<=i && Nat(i)<v.size(),"v="<<v<<" i="<<i);
    return v[i];
}


template<class X>
Vector<X> __vgetslice__(const Vector<X>& v, Int start, Int stop)
{
    if(start<0) { start+=v.size(); }
    if(stop<0) { stop+=v.size(); }
    ARIADNE_ASSERT(0<=start && start<=stop && Nat(stop)<=v.size());
    return project(v,range(start,stop));
}


template<class X>
Void __vsetitem__(Vector<X>& v, Int i, const X& x)
{
    if(i<0) { i+=v.size(); }
    ARIADNE_ASSERT(0<=i && Nat(i)<v.size());
    v[i]=x;
}


template<class X>
X __mgetitem__(const Matrix<X>& A, const boost::python::tuple& tup)
{
    Nat i=boost::python::extract<Nat>(tup[0]);
    Nat j=boost::python::extract<Nat>(tup[1]);
    return A[i][j];
}

template<class X>
Void __msetitem__(Matrix<X>& A, const boost::python::tuple& tup, const X& x)
{
    Nat i=boost::python::extract<Nat>(tup[0]);
    Nat j=boost::python::extract<Nat>(tup[1]);
    A[i][j]=x;
}

template<class X> using UniformNormType = decltype(abs(declval<X>()+declval<X>()));
template<class X> using SupremumNormType = decltype(max(abs(declval<X>()),abs(declval<X>())));
template<class X> using EuclideanNormType = decltype(sqrt(sqr(declval<X>())+sqr(declval<X>())));

template<class X> Matrix<X> transpose(const Matrix<X>& A) { return MatrixTranspose< Matrix<X> >(A); }
template<class X> UniformNormType<X> norm(const Matrix<X>& A);


template<class X>
struct from_python< Vector<X> >
{
    from_python() {
        converter::registry::push_back(&convertible,&construct,type_id< Vector<X> >());
    }

    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data)
    {
        list lst=boost::python::extract<list>(obj_ptr);
        Void* storage = ((converter::rvalue_from_python_storage< Vector<X> >*)   data)->storage.bytes;
        Vector<X> res(len(lst));
        for(Nat i=0; i!=res.size(); ++i) { res[i]=boost::python::extract<X>(lst[i]); }
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

    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static Void construct(PyObject* obj_ptr,converter::rvalue_from_python_stage1_data* data)
    {
        list rows=boost::python::extract<list>(obj_ptr);
        Matrix<X> res(len(rows),len(boost::python::extract<list>(rows[0])));
        for(Nat i=0; i!=res.row_size(); ++i) {
            list elmnts=boost::python::extract<list>(rows[i]);
            ARIADNE_ASSERT(Nat(len(elmnts))==res.column_size());
            for(Nat j=0; j!=res.column_size(); ++j) {
                res[i][j]=boost::python::extract<X>(elmnts[j]); } }
        Void* storage = ((converter::rvalue_from_python_storage< Matrix<X> >*)   data)->storage.bytes;
        new (storage) Matrix<X>(res);
        data->convertible = storage;
    }
};

OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPApproximation>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPBounds>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPValue>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloatDP>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<X>& repr) {
    return os << repr.reference(); }
template<class PR> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatValue<PR>>& repr) {
    return os << repr.reference(); }

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Vector<X> >& repr) {
    Vector<X> const& v=repr.reference();
    os << "["; for(Nat i=0; i!=v.size(); ++i) { if(i!=0) { os << ","; } os << python_representation(v[i]); } os << "]"; return os;
}

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Matrix<X> >& repr) {
    Matrix<X> const& A=repr.reference();
    os << "["; for(Nat i=0; i!=A.row_size(); ++i) { if(i!=0) { os << ","; } os << "["; for(Nat j=0; j!=A.column_size(); ++j) {
        if(j!=0) { os << ","; } os << python_representation(A[i][j]); } os << "]"; } os << "]"; return os;
}

// NOTE: g++ (version 5.4.1) does not accept decltype(auto) here
template<class V> auto __norm__(const V& v) -> decltype(norm(v)) { return norm(v); }
template<class V> auto __dot__(const V& v1, const V& v2) -> decltype(dot(v1,v2)) { return dot(v1,v2); }
template<class V1, class V2> auto __join__(const V1& v1, const V2& v2) -> decltype(join(v1,v2)) { return join(v1,v2); }
template<class X> Vector<X> __sjoin__(const X& s1, const X& s2) { return Vector<X>{s1,s2}; }
template<class M> auto __transpose__(const M& A) -> decltype(transpose(A)) { return transpose(A); }

} // namespace Ariadne


template<class X1, class X2=X1> using EqualToType = decltype(declval<X1>()==declval<X2>());
template<class X1, class X2=X1> using NotEqualToType = decltype(declval<X1>()!=declval<X2>());


template<class X>
Void export_vector_class(class_<Vector<X> >& vector_class)
{
    vector_class.def(init<Int>());
    vector_class.def(init<Int,X>());
    vector_class.def("size", &Vector<X>::size);
    vector_class.def("__len__", &Vector<X>::size);
    vector_class.def("__setitem__", &__vsetitem__<X>);
    vector_class.def("__setitem__", &__vsetitem__<X>);
    vector_class.def("__getitem__", &__vgetitem__<X>);
    vector_class.def("__getslice__", &__vgetslice__<X>);
    vector_class.def("__eq__", &__eq__<EqualToType<X,X>,Vector<X>,Vector<X> >);
    vector_class.def("__ne__", &__ne__<NotEqualToType<X,X>,Vector<X>,Vector<X> >);
    vector_class.def("__pos__", &__pos__< Vector<X>, Vector<X> >);
    vector_class.def("__neg__", &__neg__< Vector<X>, Vector<X> >);
    vector_class.def("__str__",&__cstr__< Vector<X> >);
    vector_class.def("__repr__",&__repr__< Vector<X> >);
    vector_class.def("unit",&Vector<X>::unit);
    vector_class.def("basis",&Vector<X>::basis);
    vector_class.staticmethod("unit");
    vector_class.staticmethod("basis");

    def("join", &__join__<Vector<X>,Vector<X>>);
    def("join", &__join__<Vector<X>,X>);
    def("join", &__join__<X,Vector<X>>);
    def("join", &__sjoin__<X>);

    from_python< Vector<X> >();
    to_python< Array< Vector<X> > >();
}

template<class Y, class X>
Void export_vector_conversion(class_<Vector<X> >& vector_class)
{
    //vector_class.def(init<Int,Y>());
    vector_class.def(init< Vector<Y> >());
    implicitly_convertible< Vector<Y>, Vector<X> >();
}

template<class X, class Y>
Void export_vector_arithmetic(class_<Vector<X> >& vector_class)
{
    vector_class.def("__add__",__add__< Vector<SumType<X,Y>>, Vector<X>, Vector<Y> >);
    vector_class.def("__sub__",__sub__< Vector<DifferenceType<X,Y>>, Vector<X>, Vector<Y> >);
    vector_class.def("__rmul__",__rmul__< Vector<ProductType<X,Y>>, Vector<X>, Y >);
    vector_class.def("__mul__",__mul__< Vector<QuotientType<X,Y>>, Vector<X>, Y >);
    vector_class.def("__div__",__div__< Vector<QuotientType<X,Y>>, Vector<X>, Y >);
    // Don't use boost::python style operators self+other (as below) because of
    // below) because of need to convert result expressions to Vector<R>.
    // vector_class.def(self+=Vector<Y>());
    // vector_class.def(self*=Y());
}


template<class X> Void export_vector()
{
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<X> >());
    export_vector_class<X>(vector_class);
    export_vector_arithmetic<X,X>(vector_class);
}


template<> Void export_vector<FloatDPValue>()
{
    typedef FloatDPValue X;
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<X> >());
    export_vector_class<X>(vector_class);
}

template<> Void export_vector<FloatMPValue>()
{
    typedef FloatMPValue X;
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<X> >());
    export_vector_class<X>(vector_class);
}

template<> Void export_vector<FloatDPBounds>()
{
    typedef FloatDPBounds X;
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<X> >());
    export_vector_class<X>(vector_class);
    export_vector_arithmetic<X,X>(vector_class);
    def("norm",&__norm__<Vector<X>>);
    def("dot", &__dot__<Vector<X>>);

    export_vector_conversion<FloatDPValue,FloatDPBounds>(vector_class);
}

template<> Void export_vector<FloatDPApproximation>()
{
    typedef FloatDPApproximation X;
    class_< Vector<X> > vector_class(python_name<X>("Vector"),init< Vector<FloatDPApproximation> >());
    export_vector_class<X>(vector_class);
    export_vector_arithmetic<X,X>(vector_class);
    vector_class.def("__rmul__",__rmul__< Vector<X>, Vector<X>, X >);
    vector_class.def("__mul__",__mul__< Vector<X>, Vector<X>, X >);
    vector_class.def("__div__",__div__< Vector<X>, Vector<X>, X >);

    export_vector_conversion<FloatDPBounds,FloatDPApproximation>(vector_class);
}


template<class X>
Void export_covector(class_<Covector<X> >& covector_class)
{
    covector_class.def(init<Int>());
    covector_class.def(init<Int,X>());
    covector_class.def("size", &Covector<X>::size);
    covector_class.def("__len__", &Covector<X>::size);
    covector_class.def("__setitem__", &__setitem__<Covector<X>>);
    covector_class.def("__setitem__", &__setitem__<Covector<X>>);
    covector_class.def("__getitem__", &__getitem__<Covector<X>>);
    covector_class.def("__str__",&__cstr__< Covector<X> >);

    covector_class.def("__add__",__add__< Covector<SumType<X,X>>, Covector<X>, Covector<X> >);
    covector_class.def("__sub__",__sub__< Covector<DifferenceType<X,X>>, Covector<X>, Covector<X> >);
    covector_class.def("__rmul__",__rmul__< Covector<ProductType<X,X>>, Covector<X>, X >);
    covector_class.def("__mul__",__mul__< Covector<QuotientType<X,X>>, Covector<X>, X >);
    covector_class.def("__div__",__div__< Covector<QuotientType<X,X>>, Covector<X>, X >);

    def("transpose", (Covector<X>(*)(Vector<X>const&)) &transpose);
    def("transpose", (Vector<X>(*)(Covector<X>const&)) &transpose);
}


template<class X>
Void export_matrix_class(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(init<Int,Int>());
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

    matrix_class.def("identity",(Matrix<X>(*)(SizeType)) &Matrix<X>::identity);
    matrix_class.staticmethod("identity");

    from_python< Matrix<X> >();
}


template<class Y, class X>
Void export_matrix_conversion(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(init< Matrix<Y> >());
}

template<class X, class Y>
Void export_matrix_arithmetic(class_<Matrix<X> >& matrix_class)
{
    matrix_class.def("__add__", &__add__< Matrix<SumType<X,Y>>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__sub__", &__sub__< Matrix<DifferenceType<X,Y>>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__mul__", &__mul__< Matrix<ArithmeticType<X,Y>>, Matrix<X>, Y >);
    matrix_class.def("__rmul__", &__rmul__< Matrix<ProductType<X,Y>>, Matrix<X>, Y >);
    matrix_class.def("__div__", &__div__< Matrix<QuotientType<X,Y>>, Matrix<X>, Y >);
    matrix_class.def("__mul__", &__mul__< Vector<ArithmeticType<X,Y>>, Matrix<X>, Vector<Y> >);
    matrix_class.def("__mul__", &__mul__< Matrix<ArithmeticType<X,Y>>, Matrix<X>, Matrix<Y> >);
    matrix_class.def("__rmul__", &__rmul__< Covector<ArithmeticType<X,Y>>, Matrix<X>, Covector<Y> >);
}

template<class X>
Void export_matrix_operations(class_<Matrix<X> >& matrix_class)
{
    typedef decltype(abs(declval<X>()+declval<X>())) NormType;
    def("norm",(UniformNormType<X>(*)(const Matrix<X>&)) &norm);
    def("transpose",(Matrix<X>(*)(const Matrix<X>&)) &transpose);

    def("inverse",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&)) &inverse);
    def("solve",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&,const Matrix<X>&)) &solve<X>);
    def("solve",(Vector<ArithmeticType<X>>(*)(const Matrix<X>&,const Vector<X>&)) &solve<X>);
}


template<class X> Void export_matrix()
{
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init< Matrix<X> >());
    export_matrix_class<X>(matrix_class);
    export_matrix_conversion<X,X>(matrix_class);
    export_matrix_arithmetic<X,X>(matrix_class);
    export_matrix_operations<X>(matrix_class);


}

template<> Void export_matrix<FloatDPValue>()
{
    typedef FloatDPValue X;
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init<Matrix<X>>());
    export_matrix_class<X>(matrix_class);
}

template<> Void export_matrix<FloatMPValue>()
{
    typedef FloatMPValue X;
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init<Matrix<X>>());
    export_matrix_class<X>(matrix_class);
}

template<> Void export_matrix<FloatDPBounds>()
{
    typedef FloatDPBounds X;
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init<Matrix<X>>());
    export_matrix_class<X>(matrix_class);
    export_matrix_conversion<FloatDPValue,FloatDPBounds>(matrix_class);
    export_matrix_arithmetic<X,X>(matrix_class);
    export_matrix_operations<X>(matrix_class);
    def("gs_inverse", (Matrix<X>(*)(const Matrix<X>&)) &gs_inverse);
    def("lu_inverse", (Matrix<X>(*)(const Matrix<X>&)) &lu_inverse);
    def("gs_solve", (Vector<X>(*)(const Matrix<X>&,const Vector<X>&)) &gs_solve);
    def("gs_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &gs_solve);
    def("lu_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &lu_solve);

    def("triangular_decomposition",&triangular_decomposition<X>);
    def("orthogonal_decomposition", &orthogonal_decomposition<X>);

    //implicitly_convertible< Matrix<FloatDPApproximation>, Matrix<FloatDPBounds> >();
}

template<> Void export_matrix<FloatDPApproximation>()
{
    typedef FloatDPApproximation X;
    class_< Matrix<X> > matrix_class(python_name<X>("Matrix"),init<Matrix<X>>());
    export_matrix_class<X>(matrix_class);
    export_matrix_conversion<FloatDPBounds,FloatDPApproximation>(matrix_class);
    export_matrix_arithmetic<X,X>(matrix_class);
    export_matrix_operations<X>(matrix_class);

    def("triangular_decomposition",&triangular_decomposition<X>);
    def("orthogonal_decomposition", &orthogonal_decomposition<X>);
    def("row_norms",(Vector<X>(*)(const Matrix<X>&)) &row_norms<X>);

//    to_python< Tuple<FloatMatrix,FloatMatrix,PivotMatrix> >();
}


template<class X> Void export_diagonal_matrix()
{
    class_< DiagonalMatrix<X> > diagonal_matrix_class(python_name<X>("DiagonalMatrix"),no_init);
    diagonal_matrix_class.def(init< SizeType >());
    diagonal_matrix_class.def(init< Vector<X> >());
    diagonal_matrix_class.def("__setitem__", &DiagonalMatrix<X>::set);
    diagonal_matrix_class.def("__getitem__", &DiagonalMatrix<X>::get, return_value_policy<copy_const_reference>());
    diagonal_matrix_class.def("__str__",&__cstr__< DiagonalMatrix<X> >);
    diagonal_matrix_class.def("__mul__", &__mul__< DiagonalMatrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    diagonal_matrix_class.def("__div__", &__div__< DiagonalMatrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    diagonal_matrix_class.def("__mul__", &__mul__< Vector<X>, DiagonalMatrix<X>, Vector<X> >);
    diagonal_matrix_class.def("__mul__", &__mul__< Matrix<X>, DiagonalMatrix<X>, Matrix<X> >);
    //diagonal_matrix_class.def("__rmul__", &__mul__< Covector<X>, Covector<X>, DiagonalMatrix<X> >);
    diagonal_matrix_class.def("__rmul__", &__mul__< Matrix<X>, DiagonalMatrix<X>, DiagonalMatrix<X> >);
    //diagonal_matrix_class.def("__neg__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>) operator- );
    //diagonal_matrix_class.def("__add__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator+ );
    //diagonal_matrix_class.def("__sub__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator- );
    //diagonal_matrix_class.def("__mul__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator* );
    //diagonal_matrix_class.def("__div__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator/ );
    //diagonal_matrix_class.def("__mul__", (Vector<X>(*)(const DiagonalMatrix<X>&,Vector<X>)) operator* );
    //diagonal_matrix_class.def("__mul__", (Matrix<X>(*)(const DiagonalMatrix<X>&,Matrix<X>)) operator* );
    //def("inverse", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&)) &inverse<X>);
}

Void export_pivot_matrix()
{

//    implicitly_convertible< PivotMatrix, Matrix<FloatDPValue> >();

    class_<PivotMatrix> pivot_matrix_class("PivotMatrix",no_init);
//    pivot_matrix_class.def("__str__",&__cstr__<PivotMatrix>);
//    pivot_matrix_class.def("__repr__",&__cstr__<PivotMatrix>);

}

template Void export_diagonal_matrix<FloatDPApproximation>();

template Void export_vector<Rational>();
template Void export_matrix<Rational>();


Void linear_algebra_submodule() {
    export_vector<FloatDPApproximation>();
    export_vector<FloatDPBounds>();
    export_vector<FloatDPBall>();
    export_vector<FloatDPValue>();

    export_vector<FloatMPApproximation>();
    export_vector<FloatMPBounds>();
    export_vector<FloatMPBall>();
    export_vector<FloatMPValue>();

    export_matrix<FloatDPApproximation>();
    export_matrix<FloatDPBounds>();
    export_matrix<FloatDPValue>();

    export_matrix<FloatMPApproximation>();
    export_matrix<FloatMPBounds>();
    export_matrix<FloatMPValue>();

    export_pivot_matrix();
    export_diagonal_matrix<FloatDPApproximation>();

    export_vector<Rational>();
    export_matrix<Rational>();
}
