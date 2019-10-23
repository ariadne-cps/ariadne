#pragma once

#if defined HAVE_GLPK_H
#include <bits/stdc++.h>
#include <glpk.h>
#endif

namespace Ariadne {

#if defined HAVE_GLPK_H

//  perform simplex algorithm to find minimum with only lower constriants
//    using glpk library
template<class X>
Vector<X>
lp_min(const Vector<X>& C,
       const Matrix<X>& A,
       const Vector<X>& b,
       const Vector<X>& lb,
       int& errnum);

#endif

} // namespace Ariadne

#include "glpk_interface.impl.hpp"