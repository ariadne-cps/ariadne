// Evolution algorithms and computed trajectories for FirstOrderPDE.
#ifndef ARIADNE_EVOLUTION_FIRST_ORDER_PDE_HPP
#define ARIADNE_EVOLUTION_FIRST_ORDER_PDE_HPP

#include "system/first_order_pde.hpp"
#include "numeric/float.decl.hpp"
#include "numeric/float_bounds.hpp"
#include "algebra/tensor.hpp"

namespace Ariadne {

//! \brief Compute the value of the multiaffine interpolation of the data \a us at the point \a x.
template<class X> Vector<Bounds<X>> multiaffine_interpolate(Tensor<2,Vector<Bounds<X>>> const& us, Vector<X> const& x);

//! \brief The solution of a first-order partial differential equation.
//! \relates FirstOrderPDE
template<class PR> struct FirstOrderPDESolution {
    Float<PR> h; //!< <p/>
    Float<PR> tau; //!< <p/>
    Tensor<3,Vector<FloatBounds<PR>>> uts; //!< <p/>
    FloatUpperBound<PR> error; //!< <p/>
};

//! \brief Solve the first-order partial differential equation
//! \f$ A u_{,t} + \sum_{j=1}^{n} B_j u_{,x_j} = f(u); \  u(0,x) = \phi_0(x); \  u:\R^n \to \R^m . \f$
//! \relates FirstOrderPDE
template<class PR> FirstOrderPDESolution<PR>
first_order_pde(FirstOrderPDE const& pde, EffectiveVectorMultivariateFunction const& phi0, PR pr);

} // namespace Ariadne
#endif // ARIADNE_EVOLUTION_FIRST_ORDER_PDE_HPP
