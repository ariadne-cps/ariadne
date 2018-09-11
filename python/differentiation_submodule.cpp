/***************************************************************************
 *            differentiation_submodule.cpp
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

#include "utility/typedefs.hpp"
#include "utility/array.hpp"
#include "numeric/numeric.hpp"
#include "algebra/vector.hpp"
#include "algebra/covector.hpp"
#include "algebra/differential.hpp"
#include "algebra/univariate_differential.hpp"
#include "algebra/fixed_differential.hpp"
#include "algebra/fixed_univariate_differential.hpp"
#include "algebra/expansion.inl.hpp"

using namespace Ariadne;

namespace Ariadne {

inline Nat compute_polynomial_data_size(Nat rs, Nat as, Nat d) { return rs*Ariadne::bin(d+as,as); }

template<class I, class X> Expansion<I,X> expansion_from_python(pybind11::dict pydct) {
    std::map<I,X> mp=pybind11::cast<std::map<I,X>>(pydct);
    return Expansion<I,X>(mp);
}

template<class I, class X> pybind11::dict to_python_dict(Ariadne::Expansion<I,X> const& e) {
    Nat n=e.argument_size();
    pybind11::dict res;
    pybind11::list lst;
    for(Nat i=0; i!=n; ++i) { lst.append(0); }
    I a;
    X c;
    for(typename Expansion<I,X>::ConstIterator iter=e.begin(); iter!=e.end(); ++iter) {
        a=iter->index();
        c=iter->coefficient();
        for(Nat i=0; i!=a.size(); ++i) { Int ai=a[i]; lst[i]=ai; }
        pybind11::tuple tup(lst);
        res[pybind11::object(a)]=pybind11::object(c);
    }
    return res;
}

} // namespace Ariadne


template<class DIFF>
DIFF make_dense_differential(const Nat& as, const Nat& d, const pybind11::object& obj)
{
    typedef typename DIFF::ValueType X;
    DIFF result(as,d);
    Array<X> coefficient=make_array<X>(obj);
    std::cerr<<"polynomial_data_size("<<as<<","<<d<<")="<<compute_polynomial_data_size(1u,as,d)<<"\n";
    std::cerr<<"coefficient="<<coefficient<<"\n";
    assert(coefficient.size()==compute_polynomial_data_size(1u,as,d));
    MultiIndex i(as);
    const X* ptr=coefficient.begin();
    while(i.degree()<=d) {
        //result[i]=*ptr; ++i; ++ptr;
        result->expansion().append(i,*ptr); ++i; ++ptr;
    }
    return result;
}

template<class DIFF>
DIFF
make_sparse_differential(const pybind11::object& obj,const Nat& d)
{
    typedef typename DIFF::ValueType X;
    Expansion<MultiIndex,X> expansion = pybind11::cast< Expansion<MultiIndex,X> >(obj);
    DIFF result=DIFF(expansion,d);
    return result;
}


template<class DIFF>
pybind11::list
make_differential_variables(const Nat& d, const Vector<typename DIFF::NumericType>& x)
{
    pybind11::list result;
    for(Nat i=0; i!=x.size(); ++i) {
        result.append(DIFF::variable(x.size(),d,x[i],i));
    }
    return result;
}


template<class DIFF>
Vector<DIFF>
make_differential_vector(const Nat& rs, const Nat& as, const Nat& d, const pybind11::object& obj)
{
    typedef typename DIFF::ValueType X;
    pybind11::list lst=pybind11::cast<pybind11::list>(obj);
    Array<X> coefficient = make_array<X>(obj);
    ARIADNE_ASSERT(coefficient.size()==compute_polynomial_data_size(rs,as,d));
    Vector<DIFF> result=Vector<DIFF>(rs,DIFF(as,d));
    for(SizeType i=0; i!=rs; ++i) {
        result[i]=make_sparse_differential<DIFF>(lst[i],d);
    }
    return result;
}

template<class PR> Differential<FloatBounds<PR>> make_differential_variable(SizeType n, DegreeType d, Real r, SizeType i, PR pr) {
    return Differential<FloatBounds<PR>>::variable(n,d,FloatBounds<PR>(r,pr),i);
}
template<class PR> Differential<FloatBounds<PR>> make_differential_variable(SizeType n, DegreeType d, ValidatedReal r, SizeType i, PR pr) {
    return Differential<FloatBounds<PR>>::variable(n,d,FloatBounds<PR>(r,pr),i);
}
template<class X, class Y, class PR> Vector<Differential<X>> make_variables(DegreeType d, Vector<Y> vy, PR pr) {
    return Differential<X>::variables(d,Vector<X>(vy,pr));
}



template<class C, class I, class X, EnableIf<IsSame<I,Int>> =dummy> inline
X get_item(const C& c, const I& i) {
    return c[static_cast<SizeType>(i)];
}

template<class C, class I, class X, EnableIf<Not<IsSame<I,Int>>> =dummy> inline
X get_item(const C& c, const I& i) {
    return c[i];
}

template<class C, class I, class J, class X> inline
X matrix_get_item(const C& c, const I& i, const J& j) { return c[static_cast<SizeType>(i)][j]; }

template<class C, class I, class X, EnableIf<IsSame<I,Int>> =dummy> inline
Void set_item(C& c, const I& i, const X& x) { c[static_cast<SizeType>(i)]=x; }

template<class C, class I, class X, EnableIf<Not<IsSame<I,Int>>> =dummy> inline
Void set_item(C& c, const I& i, const X& x) { c[i]=x; }

template<class C, class I, class J, class X> inline
Void matrix_set_item(C& c, const I& i, const J& j, const X& x) { c[i][j]=x; }


namespace Ariadne {

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Expansion<MultiIndex,X> >& repr);

template<class X> OutputStream& operator<<(OutputStream& os, const PythonRepresentation< Differential<X> >& repr) {
    const Differential<X>& diff=repr.reference();
    os << python_name<X>("Differential").c_str() << "(" << python_representation(diff.expansion()) << "," << diff.degree() << ")";
    //os << python_name<X>("Differential").c_str() << "(" << diff.argument_size() << "," << diff.degree() << "," << python_representation(diff.expansion()) << ")";
    return os;
}

}


template<class D> using ValueType = decltype(declval<D>().value());
template<class D> using GradientType = decltype(declval<D>().gradient());
template<class D> using HessianType = decltype(declval<D>().hessian());

template<class DIFF>
Void export_differential(pybind11::module& module, const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef typename X::GenericType Y;
    typedef DIFF D;

    static constexpr auto self = pybind11::detail::self;
    
    pybind11::class_<D> differential_class(module,name.c_str());
    differential_class.def(pybind11::init<D>());
    differential_class.def(pybind11::init(&make_sparse_differential<D>) );
    differential_class.def( pybind11::init< SizeType, DegreeType >());
    differential_class.def( pybind11::init< Expansion<MultiIndex,X>, DegreeType >());
    differential_class.def("__getitem__", &get_item<D,MultiIndex,X>);
    differential_class.def("__setitem__",&set_item<D,MultiIndex,X>);
    differential_class.def(-self);
    differential_class.def(self+self);
    differential_class.def(self-self);
    differential_class.def(self*self);
    differential_class.def(self/self);
    differential_class.def(self+X());
    differential_class.def(self-X());
    differential_class.def(self*X());
    differential_class.def(self/X());
    differential_class.def(X()+self);
    differential_class.def(X()-self);
    differential_class.def(X()*self);
    differential_class.def(self+Y());
    differential_class.def(self-Y());
    differential_class.def(self*Y());
    differential_class.def(self/Y());
    differential_class.def(Y()+self);
    differential_class.def(Y()-self);
    differential_class.def(Y()*self);
    differential_class.def(self+=self);
    differential_class.def(self-=self);
    differential_class.def(self+=X());
    differential_class.def(self-=X());
    differential_class.def(self*=X());
    differential_class.def(self/=X());
    differential_class.def("__str__", &__cstr__<D>);
    differential_class.def("__repr__", &__repr__<D>);

    differential_class.def("value", (ValueType<D>(D::*)()const)&D::value);
    differential_class.def("gradient", (GradientType<D>(D::*)()const)&D::gradient);
    differential_class.def("hessian", (HessianType<D>(D::*)()const)&D::hessian);
    differential_class.def("expansion", (Expansion<MultiIndex,X>const&(D::*)()const)&D::expansion);

    differential_class.def_static("constant",(D(*)(SizeType, DegreeType, const X&))&D::constant);
    differential_class.def_static("variable",(D(*)(SizeType, DegreeType, const X&, SizeType))&D::variable);
    differential_class.def_static("constants",(Vector<D>(*)(SizeType, DegreeType, const Vector<X>&))&D::constants);
    differential_class.def_static("variables",(Vector<D>(*)(DegreeType, const Vector<X>&))&D::variables);

    
    module.def("derivative", (D(*)(const D&, SizeType))&D::_derivative);
    module.def("antiderivative", (D(*)(const D&, SizeType))&D::_antiderivative);

    typedef D(*UFn)(D const&);
    
    module.def("neg", (D(*)(D const&)) &_neg_<D>);
    module.def("sqr", (UFn) &_sqr_<D>);
    module.def("rec",(UFn) &_rec_<D>);
    module.def("pow", (D(*)(D const&,Int const&)) &_pow_<D,Int>);
    module.def("sqrt",(UFn) &_sqrt_<D>);
    module.def("exp",(UFn) &_exp_<D>);
    module.def("log",(UFn) &_log_<D>);
    module.def("sin",(UFn) &_sin_<D>);
    module.def("cos",(UFn) &_cos_<D>);
    module.def("tan",(UFn) &_tan_<D>);
    module.def("atan",(UFn) &_atan_<D>);
}

template<class DIFF>
Void
export_differential_vector(pybind11::module& module, const String& name)
{
    typedef typename DIFF::ValueType X;
    typedef Vector<X> V;
    typedef DIFF D;
    typedef Vector<D> DV;

    pybind11::class_<DV> differential_vector_class(module,name.c_str());
    differential_vector_class.def(pybind11::init<DV>());
    differential_vector_class.def(pybind11::init(&make_differential_vector<D>));
    differential_vector_class.def(pybind11::init<Nat,Nat,Nat>());
    differential_vector_class.def("__getitem__", &matrix_get_item<DV,Int,MultiIndex,X>);
    differential_vector_class.def("__getitem__", &get_item<DV,Int,D>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,X>);
    differential_vector_class.def("__setitem__",&set_item<DV,Int,D>);
    differential_vector_class.def("__neg__",&__neg__<DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,DV>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,DV>);
    differential_vector_class.def("__add__",&__add__<DV,DV,V>);
    differential_vector_class.def("__sub__",&__sub__<DV,DV,V>);
    differential_vector_class.def("__rmul__",&__rmul__<DV,DV,D>);
    differential_vector_class.def("__mul__",&__mul__<DV,DV,D>);
    differential_vector_class.def("__div__",&__div__<DV,DV,D>);
    differential_vector_class.def("__rmul__",&__rmul__<DV,DV,X>);
    differential_vector_class.def("__mul__",&__mul__<DV,DV,X>);
    differential_vector_class.def("__div__",&__div__<DV,DV,X>);
    differential_vector_class.def("value", &DV::value);
    differential_vector_class.def("jacobian", &DV::jacobian);
    differential_vector_class.def("__str__",&__cstr__<DV>);
    //differential_vector_class.def("__repr__",&__repr__<DV>);

    module.def("compose",(D(*)(const D&,const DV&))&D::_compose);
    module.def("compose",(DV(*)(const DV&,const DV&))&DV::_compose);

    module.def("solve",(DV(*)(const DV&,const V&))&DV::_solve);
    module.def("flow",(DV(*)(const DV&,const V&))&DV::_flow);

    //module.def("lie_derivative", (DV(*)(const DV&,const DV&))&lie_derivative);
}

template Void export_differential< Differential<FloatDPApproximation> >(pybind11::module&,const String&);
template Void export_differential< Differential<FloatDPBounds> >(pybind11::module&,const String&);

template Void export_differential_vector< Differential<FloatDPApproximation> >(pybind11::module&,const String&);
template Void export_differential_vector< Differential<FloatDPBounds> >(pybind11::module&,const String&);

Void differentiation_submodule(pybind11::module& module)
{

    export_differential< Differential<FloatDPApproximation> >(module,python_name<FloatDPApproximation>("Differential"));
    export_differential< Differential<FloatDPBounds> >(module,python_name<FloatDPBounds>("Differential"));
    export_differential_vector< Differential<FloatDPApproximation> >(module,python_name<FloatDPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatDPBounds> >(module,python_name<FloatDPBounds>("DifferentialVector"));

    export_differential< Differential<FloatMPApproximation> >(module,python_name<FloatMPApproximation>("Differential"));
    export_differential< Differential<FloatMPBounds> >(module,python_name<FloatMPBounds>("Differential"));
    export_differential_vector< Differential<FloatMPApproximation> >(module,python_name<FloatMPApproximation>("DifferentialVector"));
    export_differential_vector< Differential<FloatMPBounds> >(module,python_name<FloatMPBounds>("DifferentialVector"));

    module.def("differential_variables", (Vector<Differential<FloatDPBounds>>(*)(DegreeType, Vector<Real>, DoublePrecision)) &make_variables<FloatDPBounds>);
    module.def("differential_variables", (Vector<Differential<FloatMPBounds>>(*)(DegreeType, Vector<Real>, MultiplePrecision)) &make_variables<FloatMPBounds>);


}

