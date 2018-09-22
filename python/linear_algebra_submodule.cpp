/***************************************************************************
 *            linear_algebra_submodule.cpp
 *
 *  Copyright 2008--17  Pieter Collins
 *
 ****************************************************************************/

/*
 *  This file is part of Ariadne.
 *
 *  Ariadne is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  Ariadne is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with Ariadne.  If not, see <https://www.gnu.org/licenses/>.
 */

#include "pybind11.hpp"
#include "utilities.hpp"

#include "config.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"

using namespace Ariadne;

namespace Ariadne {

template<class X>
X __vgetitem__(const Vector<X>& v, Int i)
{
    if(i<0) { i+=static_cast<Int>(v.size()); }
    ARIADNE_ASSERT_MSG(0<=i && Nat(i)<v.size(),"v="<<v<<" i="<<i);
    return v[static_cast<Nat>(i)];
}


template<class X>
Vector<X> __vgetslice__(const Vector<X>& v, Int start, Int stop)
{
    if(start<0) { start+=static_cast<Int>(v.size()); }
    if(stop<0) { stop+=static_cast<Int>(v.size()); }
    ARIADNE_ASSERT(0<=start && start<=stop && Nat(stop)<=v.size());
    return project(v,range(static_cast<Nat>(start),static_cast<Nat>(stop)));
}


template<class X>
Void __vsetitem__(Vector<X>& v, Int i, const X& x)
{
    if(i<0) { i+=static_cast<Int>(v.size()); }
    ARIADNE_ASSERT(0<=i && Nat(i)<v.size());
    v[static_cast<Nat>(i)]=x;
}


template<class X>
X __cvgetitem__(const Covector<X>& u, Int i)
{
    if(i<0) { i+=static_cast<Int>(u.size()); }
    ARIADNE_ASSERT_MSG(0<=i && Nat(i)<u.size(),"v="<<u<<" i="<<i);
    return u[static_cast<Nat>(i)];
}

template<class X>
Void __cvsetitem__(Covector<X>& u, Int i, const X& x)
{
    if(i<0) { i+=static_cast<Int>(u.size()); }
    ARIADNE_ASSERT(0<=i && Nat(i)<u.size());
    u[static_cast<Nat>(i)]=x;
}

template<class X>
X __mgetitem__(const Matrix<X>& A, const std::tuple<Nat,Nat>& ind)
{
    Nat i=std::get<0>(ind);
    Nat j=std::get<1>(ind);
    return A[i][j];
}

template<class X>
Void __msetitem__(Matrix<X>& A, const std::tuple<Nat,Nat>& ind, const X& x)
{
    Nat i=std::get<0>(ind);
    Nat j=std::get<1>(ind);
    A[i][j]=x;
}


template<class X> using UniformNormType = decltype(abs(declval<X>()+declval<X>()));
template<class X> using SupremumNormType = decltype(max(abs(declval<X>()),abs(declval<X>())));
template<class X> using EuclideanNormType = decltype(sqrt(sqr(declval<X>())+sqr(declval<X>())));

template<class X> Matrix<X> transpose(const Matrix<X>& A) { return MatrixTranspose<Matrix<X>>(A); }
template<class X> UniformNormType<X> norm(const Matrix<X>& A);


template<class X> pybind11::list vector_to_python(Vector<X> const& v) {
    pybind11::list lst;
    for(SizeType i=0; i!=v.size(); ++i) {
        lst.append(pybind11::cast(v[i]));
    }
    return lst;
}

template<class X> pybind11::list matrix_to_python(Matrix<X> const& A) {
    pybind11::list res;
    for(SizeType i=0; i!=A.row_size(); ++i) {
        pybind11::list row;
        for(SizeType j=0; j!=A.column_size(); ++j) {
            row.append(pybind11::cast(A[i][j]));
        }
        res.append(row);
    }
    return res;
}


template<class X> Vector<X> vector_from_python(pybind11::list const& lst) {
    SizeType n=lst.size();
    Vector<X> r(n);
    for(SizeType i=0; i!=n; ++i) {
        pybind11::object obj=lst[i];
        X const* val = obj.cast<X*>();
        r[i]=*val;
    }
    return r;
}

template<class X> Covector<X> covector_from_python(pybind11::object const& obj) {
    pybind11::list const& lst=static_cast<pybind11::list const&>(obj);
    SizeType n=lst.size();
    Covector<X> r(n);
    for(SizeType j=0; j!=n; ++j) {
        pybind11::object elem=lst[j];
        X const* val = elem.cast<X*>();
        r[j]=*val;
    }
    return r;
}

template<class X> Matrix<X> matrix_from_python(pybind11::list const& lst) {
    pybind11::object row0_obj=lst[0];
    pybind11::list const& row0=static_cast<pybind11::list const&>(row0_obj);
    SizeType rs=lst.size();
    SizeType cs=row0.size();
    Matrix<X> r(rs,cs);
    for(SizeType i=0; i!=rs; ++i) {
        pybind11::object row_obj=lst[i];
        pybind11::list const& row=static_cast<pybind11::list const&>(row_obj);
        for(SizeType j=0; j!=cs; ++j) {
            pybind11::object val_obj=row[j];
            X const* val = val_obj.cast<X*>();
            r[i][j]=*val;
        }
    }
    return r;
}


OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Rational>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPApproximation>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPBounds>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<FloatDPValue>& repr);
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<RawFloatDP>& repr);

template<class T> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<T>& repr) {
    return os <<repr.reference(); }

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Vector<X>>& repr) {
    Vector<X> const& v=repr.reference();
    os <<"["; for(Nat i=0; i!=v.size(); ++i) { if(i!=0) { os <<","; } os <<python_representation(v[i]); } os <<"]"; return os;
}

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Matrix<X>>& repr) {
    Matrix<X> const& A=repr.reference();
    os <<"["; for(Nat i=0; i!=A.row_size(); ++i) { if(i!=0) { os <<","; } os <<"["; for(Nat j=0; j!=A.column_size(); ++j) {
        if(j!=0) { os <<","; } os <<python_representation(A[i][j]); } os <<"]"; } os <<"]"; return os;
}

// NOTE: g++ (version 5.4.1) does not accept decltype(auto) here
template<class V> auto __norm__(const V& v) -> decltype(norm(v)) { return norm(v); }
template<class V> auto __dot__(const V& v1, const V& v2) -> decltype(dot(v1,v2)) { return dot(v1,v2); }
template<class V1, class V2> auto __join__(const V1& v1, const V2& v2) -> decltype(join(v1,v2)) { return join(v1,v2); }
template<class X> Vector<X> __sjoin__(const X& s1, const X& s2) { return Vector<X>{s1,s2}; }
template<class M> auto __transpose__(const M& A) -> decltype(transpose(A)) { return transpose(A); }

} // namespace Ariadne



template<class X>
Void define_vector_constructors(pybind11::module& module, pybind11::class_<Vector<X>>& vector_class)
{
    vector_class.def(pybind11::init<Vector<X>>());
    vector_class.def(pybind11::init<Nat>());
    vector_class.def(pybind11::init<Nat,X>());
    vector_class.def_static("unit",&Vector<X>::unit);
    vector_class.def_static("basis",&Vector<X>::basis);
    vector_class.def(pybind11::init([](pybind11::list const& lst){return Vector<X>(pybind11::cast<Array<X>>(lst));}));
    pybind11::implicitly_convertible<pybind11::list,Vector<X>>();
}

template<class Y, class X>
Void define_vector_conversion(pybind11::module& module, pybind11::class_<Vector<X>>& vector_class)
{
    vector_class.def(pybind11::init<Vector<Y>>());
    pybind11::implicitly_convertible<Vector<Y>,Vector<X>>();
}

template<class X, class Y>
Void define_mixed_vector_arithmetic(pybind11::module& module, pybind11::class_<Vector<X>>& vector_class)
{
    vector_class.def("__add__",__add__<Vector<X>,Vector<Y> , Return<Vector<SumType<X,Y>>>>);
    vector_class.def("__sub__",__sub__<Vector<X>,Vector<Y> , Return<Vector<DifferenceType<X,Y>>> >);
    vector_class.def("__rmul__",__rmul__<Vector<X>,Y , Return<Vector<ProductType<Y,X>>> >);
    vector_class.def("__mul__",__mul__<Vector<X>,Y , Return<Vector<ProductType<X,Y>>> >);
    vector_class.def("__div__",__div__<Vector<X>,Y , Return<Vector<QuotientType<X,Y>>> >);
    module.def("dot",  &_dot_<Vector<X>,Vector<Y>>);
    // Don't use operators self+other (as below) because of
    // below) because of need to convert result expressions to Vector<R>.
    // vector_class.def(self+=Vector<Y>());
    // vector_class.def(self*=Y());
}

template<class X>
Void define_vector(pybind11::module& module, pybind11::class_<Vector<X>>& vector_class)
{
    define_vector_constructors(module,vector_class);
    define_vector_concept(module,vector_class);

//    module.def("join", &__join__<Vector<X>,Vector<X>>);
//    module.def("join", &__join__<Vector<X>,X>);
//    module.def("join", &__join__<X,Vector<X>>);
//    module.def("join", &__sjoin__<X>);

//    vector_class.def(pybind11::init(&vector_from_python<X>));
//    pybind11::implicitly_convertible<pybind11::list, Vector<X>>();
}


template<class F> Void define_vector(pybind11::module& module, pybind11::class_<Vector<Value<F>>>& vector_class) {
    define_vector_constructors(module, vector_class);
    define_vector_operations(module, vector_class);
}

template<class F> Void define_vector(pybind11::module& module, pybind11::class_<Vector<Bounds<F>>>& vector_class) {
    using X=Bounds<F>;
    using PR=PrecisionType<X>;
    define_vector_constructors(module, vector_class);
    define_vector_concept(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);

//    vector_class.def(pybind11::init<Vector<ValidatedNumber>,PR>());
//    vector_class.def(pybind11::init<Array<ValidatedNumber>,PR>());
    vector_class.def(pybind11::init<Vector<Rational>,PR>());
    vector_class.def(pybind11::init<Vector<Value<F>>>());
    pybind11::implicitly_convertible<Vector<Value<F>>,Vector<Bounds<F>>>();
}

template<class F> Void define_vector(pybind11::module& module, pybind11::class_<Vector<Approximation<F>>>& vector_class)
{
    using X=Approximation<F>;
    using PR=PrecisionType<X>;
    define_vector_constructors(module, vector_class);
    define_vector_concept(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);

//    vector_class.def(pybind11::init<Vector<ApproximateNumber>,PR>());
//    vector_class.def(pybind11::init<Array<ApproximateNumber>,PR>());
    vector_class.def(pybind11::init<Vector<Rational>,PR>());
    vector_class.def(pybind11::init<Vector<Bounds<F>>>());
    pybind11::implicitly_convertible<Vector<Bounds<F>>,Vector<Approximation<F>>>();
}


template<class X> Void export_vector(pybind11::module& module)
{
    pybind11::class_<Vector<X>> vector_class(module,(class_name<X>()+"Vector").c_str());
    define_vector(module, vector_class);
}


template<class X>
Void define_covector(pybind11::module& module, pybind11::class_<Covector<X>>& covector_class)
{
    covector_class.def(pybind11::init<Covector<X>>());
    covector_class.def(pybind11::init<Int>());
    covector_class.def(pybind11::init<Int,X>());
    covector_class.def("size", &Covector<X>::size);
    covector_class.def("__len__", &Covector<X>::size);
    covector_class.def("__setitem__", &__setitem__<Covector<X>,Int,X>);
    covector_class.def("__getitem__", &__getitem__<Covector<X>,Int>);
    covector_class.def("__str__",&__cstr__<Covector<X>>);

    covector_class.def("__add__",__add__<Covector<X>,Covector<X> , Return<Covector<SumType<X,X>>> >);
    covector_class.def("__sub__",__sub__<Covector<X>,Covector<X> , Return<Covector<DifferenceType<X,X>>> >);
    covector_class.def("__rmul__",__rmul__<Covector<X>,X , Return<Covector<ProductType<X,X>>> >);
    covector_class.def("__mul__",__mul__<Covector<X>,X , Return<Covector<QuotientType<X,X>>> >);
    covector_class.def("__div__",__div__<Covector<X>,X , Return<Covector<QuotientType<X,X>>> >);

    module.def("transpose", (Covector<X>const&(*)(Vector<X>const&)) &transpose);
    module.def("transpose", (Vector<X>const&(*)(Covector<X>const&)) &transpose);
}

template<class X>
Void export_covector(pybind11::module& module)
{
    pybind11::class_<Covector<X>> covector_class(module,python_name<X>("Covector").c_str());
    define_covector(module,covector_class);
}

template<class X> class Tag { };

template<class X>
Void define_matrix_class(pybind11::module& module, pybind11::class_<Matrix<X>>& matrix_class)
{
    matrix_class.def(pybind11::init<Matrix<X>>());
    if constexpr (IsDefaultConstructible<X>::value) {
        matrix_class.def(pybind11::init<Int,Int>()); }
    matrix_class.def(pybind11::init<Int,Int,X>());
    matrix_class.def("rows", &Matrix<X>::row_size);
    matrix_class.def("columns", &Matrix<X>::column_size);
    matrix_class.def("row_size", &Matrix<X>::row_size);
    matrix_class.def("column_size", &Matrix<X>::column_size);
    matrix_class.def("__setitem__", &__msetitem__<X>);
    matrix_class.def("__getitem__", &__mgetitem__<X>);
    matrix_class.def("__str__",&__cstr__<Matrix<X>>);
    matrix_class.def("__repr__",&__repr__<Matrix<X>>);

    matrix_class.def_static("identity",(Matrix<X>(*)(SizeType)) &Matrix<X>::identity);

    matrix_class.def(pybind11::init([](pybind11::list const& lst){return matrix_from_python<X>(lst);}));
    pybind11::implicitly_convertible<pybind11::list,Matrix<X>>();
}


template<class Y, class X>
Void define_matrix_conversion(pybind11::module& module, pybind11::class_<Matrix<X>>& matrix_class)
{
    matrix_class.def(pybind11::init<Matrix<Y>>());
}

template<class X, class Y>
Void define_matrix_arithmetic(pybind11::module& module, pybind11::class_<Matrix<X>>& matrix_class)
{
    matrix_class.def("__pos__", &__pos__<Matrix<X> , Return<Matrix<X>> >);
    matrix_class.def("__neg__", &__neg__<Matrix<X> , Return<Matrix<X>> >);
    matrix_class.def("__add__", &__add__<Matrix<X>,Matrix<Y> , Return<Matrix<SumType<X,Y>>> >);
    matrix_class.def("__sub__", &__sub__<Matrix<X>,Matrix<Y> , Return<Matrix<DifferenceType<X,Y>>> >);
    matrix_class.def("__mul__", &__mul__<Matrix<X>,Y , Return<Matrix<ArithmeticType<X,Y>>> >);
    matrix_class.def("__rmul__", &__rmul__<Matrix<X>,Y , Return<Matrix<ProductType<X,Y>>> >);
    matrix_class.def("__div__", &__div__<Matrix<X>,Y , Return<Matrix<QuotientType<X,Y>>> >);
    matrix_class.def("__mul__", &__mul__<Matrix<X>,Vector<Y> , Return<Vector<ArithmeticType<X,Y>>> >);
    matrix_class.def("__mul__", &__mul__<Matrix<X>,Matrix<Y> , Return<Matrix<ArithmeticType<X,Y>>> >);
    matrix_class.def("__rmul__", &__rmul__<Matrix<X>,Covector<Y> , Return<Covector<ArithmeticType<X,Y>>> >);
}

template<class X>
Void define_matrix_operations(pybind11::module& module, pybind11::class_<Matrix<X>>& matrix_class)
{
    module.def("norm",(UniformNormType<X>(*)(const Matrix<X>&)) &norm);
    module.def("transpose",(Matrix<X>(*)(const Matrix<X>&)) &transpose);

    module.def("inverse",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&)) &inverse);
    module.def("solve",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&,const Matrix<X>&)) &solve<X>);
    module.def("solve",(Vector<ArithmeticType<X>>(*)(const Matrix<X>&,const Vector<X>&)) &solve<X>);
}



template<class F> Void define_matrix(pybind11::module& module, pybind11::class_<Matrix<Value<F>>>& matrix_class)
{
    using X=Value<F>;
    define_matrix_class<X>(module,matrix_class);
}

template<class F> Void define_matrix(pybind11::module& module, pybind11::class_<Matrix<Bounds<F>>>& matrix_class)
{
    using X=Bounds<F>;
    using PR=PrecisionType<X>;
    define_matrix_class<X>(module,matrix_class);
    matrix_class.def(pybind11::init<Matrix<Rational>, PR>());
    define_matrix_conversion<FloatValue<PR>,FloatBounds<PR>>(module,matrix_class);
    define_matrix_arithmetic<X,X>(module,matrix_class);
    define_matrix_operations<X>(module,matrix_class);
    module.def("gs_inverse", (Matrix<X>(*)(const Matrix<X>&)) &gs_inverse);
    module.def("lu_inverse", (Matrix<X>(*)(const Matrix<X>&)) &lu_inverse);
    module.def("gs_solve", (Vector<X>(*)(const Matrix<X>&,const Vector<X>&)) &gs_solve);
    module.def("gs_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &gs_solve);
    module.def("lu_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &lu_solve);

    module.def("triangular_decomposition",&triangular_decomposition<X>);
    module.def("orthogonal_decomposition", &orthogonal_decomposition<X>);

    //pybind11::implicitly_convertible<Matrix<FloatDPApproximation>, Matrix<FloatDPBounds>>();
}

template<class F> Void define_matrix(pybind11::module& module, pybind11::class_<Matrix<Approximation<F>>>& matrix_class)
{
    using X=Approximation<F>;
    using PR=PrecisionType<X>;
    define_matrix_class<X>(module,matrix_class);
    matrix_class.def(pybind11::init<Matrix<Rational>, PR>());
    define_matrix_conversion<FloatBounds<PR>,FloatApproximation<PR>>(module,matrix_class);
    define_matrix_arithmetic<X,X>(module,matrix_class);
    define_matrix_operations<X>(module,matrix_class);

    module.def("triangular_decomposition",&triangular_decomposition<X>);
    module.def("orthogonal_decomposition", &orthogonal_decomposition<X>);
    module.def("row_norms",(Vector<X>(*)(const Matrix<X>&)) &row_norms<X>);

//    to_python<Tuple<FloatMatrix,FloatMatrix,PivotMatrix>>();
}

template<class X> Void define_matrix(pybind11::module& module, pybind11::class_<Matrix<X>>& matrix_class)
{
    matrix_class.def(pybind11::init<Matrix<X>>());
    define_matrix_class<X>(module,matrix_class);
    define_matrix_conversion<X,X>(module,matrix_class);
    define_matrix_arithmetic<X,X>(module,matrix_class);
    define_matrix_operations<X>(module,matrix_class);
}


template<class X> Void export_matrix(pybind11::module& module)
{
    pybind11::class_<Matrix<X>> matrix_class(module,python_name<X>("Matrix").c_str());
    define_matrix(module,matrix_class);
}

template<class X> Void export_diagonal_matrix(pybind11::module& module)
{
    pybind11::class_<DiagonalMatrix<X>> diagonal_matrix_class(module,python_name<X>("DiagonalMatrix").c_str());
    diagonal_matrix_class.def(pybind11::init<SizeType >());
    diagonal_matrix_class.def(pybind11::init<Vector<X>>());
    diagonal_matrix_class.def("__setitem__", &DiagonalMatrix<X>::set);
    diagonal_matrix_class.def("__getitem__", &DiagonalMatrix<X>::get);
    diagonal_matrix_class.def("__str__",&__cstr__<DiagonalMatrix<X>>);
    diagonal_matrix_class.def("__mul__", &__mul__<DiagonalMatrix<X>,DiagonalMatrix<X> , Return<DiagonalMatrix<X>> >);
    diagonal_matrix_class.def("__div__", &__div__<DiagonalMatrix<X>,DiagonalMatrix<X> , Return<DiagonalMatrix<X>> >);
    diagonal_matrix_class.def("__mul__", &__mul__<DiagonalMatrix<X>,Vector<X> , Return<Vector<X>> >);
    diagonal_matrix_class.def("__mul__", &__mul__<DiagonalMatrix<X>,Matrix<X> , Return<Matrix<X>> >);
    //diagonal_matrix_class.def("__rmul__", &__mul__<Covector<X>,DiagonalMatrix<X> , Return<Covector<X>> >);
    diagonal_matrix_class.def("__rmul__", &__mul__<DiagonalMatrix<X>,DiagonalMatrix<X> , Return<Matrix<X>> >);
    //diagonal_matrix_class.def("__neg__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>) operator- );
    //diagonal_matrix_class.def("__add__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator+ );
    //diagonal_matrix_class.def("__sub__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator- );
    //diagonal_matrix_class.def("__mul__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator* );
    //diagonal_matrix_class.def("__div__", (DiagonalMatrix<X>(*)(DiagonalMatrix<X>,const DiagonalMatrix<X>&)) operator/ );
    //diagonal_matrix_class.def("__mul__", (Vector<X>(*)(const DiagonalMatrix<X>&,Vector<X>)) operator* );
    //diagonal_matrix_class.def("__mul__", (Matrix<X>(*)(const DiagonalMatrix<X>&,Matrix<X>)) operator* );
    //def("inverse", (DiagonalMatrix<X>(*)(const DiagonalMatrix<X>&)) &inverse<X>);
}

Void export_pivot_matrix(pybind11::module& module)
{

//    pybind11::implicitly_convertible<PivotMatrix, Matrix<FloatDPValue>>();

    pybind11::class_<PivotMatrix> pivot_matrix_class(module,"PivotMatrix");
//    pivot_matrix_class.def("__str__",&__cstr__<PivotMatrix>);
//    pivot_matrix_class.def("__repr__",&__cstr__<PivotMatrix>);

}


Void linear_algebra_submodule(pybind11::module& module) {
    export_vector<FloatDPApproximation>(module);
    export_vector<FloatDPBounds>(module);
    export_vector<FloatDPBall>(module);
    export_vector<FloatDPValue>(module);

    export_vector<FloatMPApproximation>(module);
    export_vector<FloatMPBounds>(module);
    export_vector<FloatMPBall>(module);
    export_vector<FloatMPValue>(module);

    export_covector<FloatDPApproximation>(module);
    export_covector<FloatDPBounds>(module);
    export_covector<FloatMPApproximation>(module);
    export_covector<FloatMPBounds>(module);

    export_matrix<FloatDPApproximation>(module);
    export_matrix<FloatDPBounds>(module);
    export_matrix<FloatDPValue>(module);

    export_matrix<FloatMPApproximation>(module);
    export_matrix<FloatMPBounds>(module);
    export_matrix<FloatMPValue>(module);

    export_pivot_matrix(module);
    export_diagonal_matrix<FloatDPApproximation>(module);

    export_vector<Real>(module);

    export_vector<Rational>(module);
    export_covector<Rational>(module);
    export_matrix<Rational>(module);
}
