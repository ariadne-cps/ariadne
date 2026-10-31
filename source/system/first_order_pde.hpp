// Mathematical model only. Solvers and solution histories live in evolution/.
#ifndef ARIADNE_SYSTEM_FIRST_ORDER_PDE_HPP
#define ARIADNE_SYSTEM_FIRST_ORDER_PDE_HPP

#include "numeric/float.decl.hpp"
#include "numeric/real.hpp"
#include "numeric/rational.hpp"
#include "algebra/matrix.hpp"
#include "algebra/tensor.hpp"
#include "function/function.decl.hpp"

namespace Ariadne {

//! \brief The first-order partial differential equation
//! \f$ A u_{,t} + \sum_{j=1}^{n} B_j u_{,x_j} = f(u). \f$
struct FirstOrderPDE {
    Matrix<Real> A; //!< <p/>
    Array<Matrix<Real>> Bs; //!< <p/>
    Array<DiagonalMatrix<Real>> const& Ds; //!< <p/>
    Array<Matrix<Real>> const& Ts; //!< <p/>
    EffectiveVectorMultivariateFunction const& f; //!< <p/>
};

// A, Bs are owned by this value. Ds, Ts and f are borrowed references:
// the caller must keep them alive for every solver invocation.

} // namespace Ariadne
#endif // ARIADNE_SYSTEM_FIRST_ORDER_PDE_HPP
