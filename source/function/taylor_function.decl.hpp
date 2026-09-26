/***************************************************************************
 *            function/taylor_function.decl.hpp
 *
 *  Forward declarations for Taylor function factories.
 ****************************************************************************/

#ifndef ARIADNE_TAYLOR_FUNCTION_DECL_HPP
#define ARIADNE_TAYLOR_FUNCTION_DECL_HPP

#include "utility/declarations.hpp"

namespace Ariadne {

template<class F> class Sweeper;
template<class P> class FunctionPatchFactoryInterface;

FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory();
FunctionPatchFactoryInterface<ValidatedTag>* make_taylor_function_patch_factory(Sweeper<FloatDP> const& sweeper);

} // namespace Ariadne

#endif // ARIADNE_TAYLOR_FUNCTION_DECL_HPP
