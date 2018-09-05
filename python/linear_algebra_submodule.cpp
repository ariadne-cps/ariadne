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
#include <pybind11/stl.h>

#include "utilities.hpp"

#include "config.hpp"

#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/matrix.hpp"
#include "algebra/diagonal_matrix.hpp"

using namespace Ariadne;


namespace Ariadne {

template<class X1, class X2> ArithmeticType<X1,X2> dot(const Vector<X1>& v1, const Vector<X2>& v2);

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

template<class X> Matrix<X> transpose(const Matrix<X>& A) { return MatrixTranspose< Matrix<X> >(A); }
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

    
/*

template<class X>
struct from_python< Vector<X> >
{
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Vector<X> >());
    }

    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        boost::python::list lst=boost::python::extract<boost::python::list>(obj_ptr);
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Vector<X> >*)   data)->storage.bytes;
        Vector<X> res(static_cast<SizeType>(len(lst)));
        for(Nat i=0; i!=res.size(); ++i) { res[i]=boost::python::extract<X>(lst[i]); }
        new (storage) Vector<X>(res);
        data->convertible = storage;
    }
};

template<class X>
struct from_python< Matrix<X> >
{
    from_python() {
        boost::python::converter::registry::push_back(&convertible,&construct,boost::python::type_id< Matrix<X> >());
    }

    static Void* convertible(PyObject* obj_ptr) {
        if (!PyList_Check(obj_ptr)) { return 0; }
        return obj_ptr;
    }

    static Void construct(PyObject* obj_ptr,boost::python::converter::rvalue_from_python_stage1_data* data)
    {
        boost::python::list rows=boost::python::extract<boost::python::list>(obj_ptr);
        Matrix<X> res(static_cast<SizeType>(len(rows)),static_cast<SizeType>(len(boost::python::extract<boost::python::list>(rows[0]))));
        for(Nat i=0; i!=res.row_size(); ++i) {
            boost::python::list elmnts=boost::python::extract<boost::python::list>(rows[i]);
            ARIADNE_ASSERT(Nat(len(elmnts))==res.column_size());
            for(Nat j=0; j!=res.column_size(); ++j) {
                res[i][j]=boost::python::extract<X>(elmnts[j]); } }
        Void* storage = ((boost::python::converter::rvalue_from_python_storage< Matrix<X> >*)   data)->storage.bytes;
        new (storage) Matrix<X>(res);
        data->convertible = storage;
    }
};

*/


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
OutputStream& operator<<(OutputStream& os, const PythonRepresentation<Real>& repr) {
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
Void export_vector_class(pybind11::module& module, pybind11::class_<Vector<X> >& vector_class)
{
    vector_class.def(pybind11::init<Int>());
    vector_class.def(pybind11::init<Int,X>());
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
    vector_class.def_static("unit",&Vector<X>::unit);
    vector_class.def_static("basis",&Vector<X>::basis);

    module.def("join", &__join__<Vector<X>,Vector<X>>);
    module.def("join", &__join__<Vector<X>,X>);
    module.def("join", &__join__<X,Vector<X>>);
    module.def("join", &__sjoin__<X>);

//    from_python< Vector<X> >();
//    to_python< Array< Vector<X> > >();
}

template<class Y, class X>
Void export_vector_conversion(pybind11::module& module, pybind11::class_<Vector<X> >& vector_class)
{
    //vector_class.def(pybind11::init<Int,Y>());
    vector_class.def(pybind11::init< Vector<Y> >());
    pybind11::implicitly_convertible< Vector<Y>, Vector<X> >();
}

template<class X, class Y=X>
Void export_vector_arithmetic(pybind11::module& module, pybind11::class_<Vector<X> >& vector_class)
{
    vector_class.def("__add__",__add__< Vector<SumType<X,Y>>, Vector<X>, Vector<Y> >);
    vector_class.def("__sub__",__sub__< Vector<DifferenceType<X,Y>>, Vector<X>, Vector<Y> >);
    vector_class.def("__rmul__",__rmul__< Vector<ProductType<X,Y>>, Vector<X>, Y >);
    vector_class.def("__mul__",__mul__< Vector<ProductType<X,Y>>, Vector<X>, Y >);
    vector_class.def("__div__",__div__< Vector<QuotientType<X,Y>>, Vector<X>, Y >);
    module.def("dot", (ArithmeticType<X,Y>(*)(Vector<X> const&, Vector<Y> const&)) &dot);
    // Don't use operators self+other (as below) because of
    // below) because of need to convert result expressions to Vector<R>.
    // vector_class.def(self+=Vector<Y>());
    // vector_class.def(self*=Y());
}


template<class X> Void export_vector(pybind11::module& module)
{
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<X> >());
    export_vector_class<X>(module, vector_class);
    export_vector_arithmetic<X,X>(module, vector_class);
}

template<> Void export_vector<FloatDPValue>(pybind11::module& module)
{
    typedef FloatDPValue X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<X> >());
    export_vector_class<X>(module, vector_class);
}

template<> Void export_vector<FloatMPValue>(pybind11::module& module)
{
    typedef FloatMPValue X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<X> >());
    export_vector_class<X>(module, vector_class);
}

template<> Void export_vector<FloatDPBounds>(pybind11::module& module)
{
    typedef FloatDPBounds X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<X> >());
    export_vector_conversion<FloatDPValue,FloatDPBounds>(module, vector_class);
    vector_class.def(pybind11::init<Vector<ValidatedNumber>,DoublePrecision>());
    export_vector_class<X>(module, vector_class);
    export_vector_arithmetic<X,X>(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);
}

template<> Void export_vector<FloatMPBounds>(pybind11::module& module)
{
    typedef FloatMPBounds X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<X> >());
    vector_class.def(pybind11::init<Vector<ValidatedNumber>,MultiplePrecision>());
    export_vector_conversion<FloatMPValue,FloatMPBounds>(module, vector_class);
    export_vector_class<X>(module, vector_class);
    export_vector_arithmetic<X,X>(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);
}

template<> Void export_vector<FloatDPApproximation>(pybind11::module& module)
{
    typedef FloatDPApproximation X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<FloatDPApproximation> >());
    vector_class.def(pybind11::init<Vector<ApproximateNumber>, DoublePrecision>());
    export_vector_class<X>(module, vector_class);
    export_vector_arithmetic<X,X>(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);
    //export_vector_conversion<FloatDPBounds,FloatDPApproximation>(module, vector_class);
    vector_class.def(pybind11::init< Vector<FloatDPBounds> >());
}

template<> Void export_vector<FloatMPApproximation>(pybind11::module& module)
{
    typedef FloatMPApproximation X;
    pybind11::class_< Vector<X> > vector_class(module,python_name<X>("Vector").c_str());
    vector_class.def(pybind11::init< Vector<FloatMPApproximation> >());
    vector_class.def(pybind11::init<Vector<ApproximateNumber>,MultiplePrecision>());
    export_vector_class<X>(module, vector_class);
    export_vector_arithmetic<X,X>(module, vector_class);
    module.def("norm",&__norm__<Vector<X>>);
    //export_vector_conversion<FloatMPBounds,FloatMPApproximation>(module, vector_class);
    vector_class.def(pybind11::init< Vector<FloatMPBounds> >());
}


template<class X>
Void export_covector(pybind11::module& module, pybind11::class_<Covector<X>>& covector_class)
{
    covector_class.def(pybind11::init<Int>());
    covector_class.def(pybind11::init<Int,X>());
    covector_class.def("size", &Covector<X>::size);
    covector_class.def("__len__", &Covector<X>::size);
//    covector_class.def("__setitem__", &__setitem__<Covector<X>>);
//    covector_class.def("__setitem__", &__setitem__<Covector<X>>);
//      covector_class.def("__getitem__", &__getitem__<Covector<X>>);
    covector_class.def("__str__",&__cstr__< Covector<X> >);

    covector_class.def("__add__",__add__< Covector<SumType<X,X>>, Covector<X>, Covector<X> >);
    covector_class.def("__sub__",__sub__< Covector<DifferenceType<X,X>>, Covector<X>, Covector<X> >);
    covector_class.def("__rmul__",__rmul__< Covector<ProductType<X,X>>, Covector<X>, X >);
    covector_class.def("__mul__",__mul__< Covector<QuotientType<X,X>>, Covector<X>, X >);
    covector_class.def("__div__",__div__< Covector<QuotientType<X,X>>, Covector<X>, X >);

    module.def("transpose", (Covector<X>const&(*)(Vector<X>const&)) &transpose);
    module.def("transpose", (Vector<X>const&(*)(Covector<X>const&)) &transpose);
}

template<class X>
Void export_covector(pybind11::module& module)
{
    pybind11::class_<Covector<X>> covector_class(module,python_name<X>("Covector").c_str());
    export_covector(module,covector_class);
}


template<class X>
Void export_matrix_class(pybind11::module& module, pybind11::class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(pybind11::init<Int,Int>());
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

    matrix_class.def_static("identity",(Matrix<X>(*)(SizeType)) &Matrix<X>::identity);

//    from_python< Matrix<X> >();
}


template<class Y, class X>
Void export_matrix_conversion(pybind11::module& module, pybind11::class_<Matrix<X> >& matrix_class)
{
    matrix_class.def(pybind11::init< Matrix<Y> >());
}

template<class X, class Y>
Void export_matrix_arithmetic(pybind11::module& module, pybind11::class_<Matrix<X> >& matrix_class)
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
Void export_matrix_operations(pybind11::module& module, pybind11::class_<Matrix<X> >& matrix_class)
{
    module.def("norm",(UniformNormType<X>(*)(const Matrix<X>&)) &norm);
    module.def("transpose",(Matrix<X>(*)(const Matrix<X>&)) &transpose);

    module.def("inverse",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&)) &inverse);
    module.def("solve",(Matrix<ArithmeticType<X>>(*)(const Matrix<X>&,const Matrix<X>&)) &solve<X>);
    module.def("solve",(Vector<ArithmeticType<X>>(*)(const Matrix<X>&,const Vector<X>&)) &solve<X>);
}


template<class X> Void export_matrix(pybind11::module& module)
{
    pybind11::class_< Matrix<X> > matrix_class(module,python_name<X>("Matrix").c_str());
    matrix_class.def(pybind11::init< Matrix<X> >());
    export_matrix_class<X>(module,matrix_class);
    export_matrix_conversion<X,X>(module,matrix_class);
    export_matrix_arithmetic<X,X>(module,matrix_class);
    export_matrix_operations<X>(module,matrix_class);


}

template<> Void export_matrix<FloatDPValue>(pybind11::module& module)
{
    typedef FloatDPValue X;
    pybind11::class_< Matrix<X> > matrix_class(module,python_name<X>("Matrix").c_str());
    matrix_class.def(pybind11::init<Matrix<X>>());
    export_matrix_class<X>(module,matrix_class);
}

template<> Void export_matrix<FloatMPValue>(pybind11::module& module)
{
    typedef FloatMPValue X;
    pybind11::class_< Matrix<X> > matrix_class(module,python_name<X>("Matrix").c_str());
    matrix_class.def(pybind11::init<Matrix<X>>());
    export_matrix_class<X>(module,matrix_class);
}

template<class PR> Void export_float_bounds_matrix(pybind11::module& module)
{
    typedef FloatBounds<PR> X;
    pybind11::class_< Matrix<X> > matrix_class(module,python_name<X>("Matrix").c_str());
    matrix_class.def(pybind11::init<Matrix<X>>());
    matrix_class.def(pybind11::init<Matrix<Rational>, PR>());
    export_matrix_class<X>(module,matrix_class);
    export_matrix_conversion<FloatValue<PR>,FloatBounds<PR>>(module,matrix_class);
    export_matrix_arithmetic<X,X>(module,matrix_class);
    export_matrix_operations<X>(module,matrix_class);
    module.def("gs_inverse", (Matrix<X>(*)(const Matrix<X>&)) &gs_inverse);
    module.def("lu_inverse", (Matrix<X>(*)(const Matrix<X>&)) &lu_inverse);
    module.def("gs_solve", (Vector<X>(*)(const Matrix<X>&,const Vector<X>&)) &gs_solve);
    module.def("gs_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &gs_solve);
    module.def("lu_solve", (Matrix<X>(*)(const Matrix<X>&,const Matrix<X>&)) &lu_solve);

    module.def("triangular_decomposition",&triangular_decomposition<X>);
    module.def("orthogonal_decomposition", &orthogonal_decomposition<X>);

    //pybind11::implicitly_convertible< Matrix<FloatDPApproximation>, Matrix<FloatDPBounds> >();
}

template<> Void export_matrix<FloatBounds<DP>>(pybind11::module& module) { export_float_bounds_matrix<DP>(module); }
template<> Void export_matrix<FloatBounds<MP>>(pybind11::module& module) { export_float_bounds_matrix<MP>(module); }

template<class PR> Void export_float_approximation_matrix(pybind11::module& module)
{
    typedef FloatApproximation<PR> X;
    pybind11::class_< Matrix<X> > matrix_class(module,python_name<X>("Matrix").c_str());
    matrix_class.def(pybind11::init<Matrix<X>>());
    export_matrix_class<X>(module,matrix_class);
    matrix_class.def(pybind11::init<Matrix<Rational>, PR>());
    export_matrix_conversion<FloatBounds<PR>,FloatApproximation<PR>>(module,matrix_class);
    export_matrix_arithmetic<X,X>(module,matrix_class);
    export_matrix_operations<X>(module,matrix_class);

    module.def("triangular_decomposition",&triangular_decomposition<X>);
    module.def("orthogonal_decomposition", &orthogonal_decomposition<X>);
    module.def("row_norms",(Vector<X>(*)(const Matrix<X>&)) &row_norms<X>);

//    to_python< Tuple<FloatMatrix,FloatMatrix,PivotMatrix> >();
}

template<> Void export_matrix<FloatApproximation<DP>>(pybind11::module& module) { export_float_approximation_matrix<DP>(module); }
template<> Void export_matrix<FloatApproximation<MP>>(pybind11::module& module) { export_float_approximation_matrix<MP>(module); }


template<class X> Void export_diagonal_matrix(pybind11::module& module)
{
    pybind11::class_< DiagonalMatrix<X> > diagonal_matrix_class(module,python_name<X>("DiagonalMatrix").c_str());
    diagonal_matrix_class.def(pybind11::init< SizeType >());
    diagonal_matrix_class.def(pybind11::init< Vector<X> >());
    diagonal_matrix_class.def("__setitem__", &DiagonalMatrix<X>::set);
    diagonal_matrix_class.def("__getitem__", &DiagonalMatrix<X>::get);
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

Void export_pivot_matrix(pybind11::module& module)
{

//    pybind11::implicitly_convertible< PivotMatrix, Matrix<FloatDPValue> >();

    pybind11::class_<PivotMatrix> pivot_matrix_class(module,"PivotMatrix");
//    pivot_matrix_class.def("__str__",&__cstr__<PivotMatrix>);
//    pivot_matrix_class.def("__repr__",&__cstr__<PivotMatrix>);

}

template Void export_diagonal_matrix<FloatDPApproximation>(pybind11::module& module);

template Void export_vector<Rational>(pybind11::module& module);
template Void export_matrix<Rational>(pybind11::module& module);


Void linear_algebra_submodule(pybind11::module& module) {
//    from_python<Vector<ApproximateNumber>>();
//    from_python<Vector<ValidatedNumber>>();

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
    export_matrix<Rational>(module);

//    pybind11::implicitly_convertible<Vector<Real>,Vector<ApproximateNumber>>();
//    pybind11::implicitly_convertible<Vector<Real>,Vector<ValidatedNumber>>();
}
