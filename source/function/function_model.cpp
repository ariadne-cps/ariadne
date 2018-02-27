/***************************************************************************
 *            function_model.cpp
 *
 *  Copyright 2011--17  Pieter Collins
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

#include "function/function_model.hpp"
#include "function/function_model.tpl.hpp"
#include "function/taylor_model.hpp"

#include "algebra/algebra.hpp"

#include "function/formula.hpp"

namespace Ariadne {

template class FunctionModelFactoryInterface<ValidatedTag,DoublePrecision,DoublePrecision>;
template class FunctionModelFactoryInterface<ValidatedTag,MultiplePrecision,MultiplePrecision>;

template class ScalarFunctionModel<ValidatedTag,IntervalDomainType,DoublePrecision>;
template class ScalarFunctionModel<ValidatedTag,IntervalDomainType,MultiplePrecision>;
template class VectorFunctionModel<ValidatedTag,IntervalDomainType,DoublePrecision>;
template class VectorFunctionModel<ValidatedTag,IntervalDomainType,MultiplePrecision>;

template class ScalarFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;
template class ScalarFunctionModel<ValidatedTag,BoxDomainType,MultiplePrecision>;
template class VectorFunctionModel<ValidatedTag,BoxDomainType,DoublePrecision>;
template class VectorFunctionModel<ValidatedTag,BoxDomainType,MultiplePrecision>;

} // namespace Ariadne
