/***************************************************************************
 *            polynomial.cc
 *
 *  Copyright 2008-15  Pieter Collins
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
 *  GNU Library General Public License for more detai1ls.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.
 */

#include "numeric/numeric.h"
#include "config.h"

#include "geometry/interval.h"
#include "function/polynomial.h"
#include "function/polynomial.tcc"

namespace Ariadne {

template class Polynomial<Float>;
template class Polynomial<ApproximateFloat>;
template class Polynomial<ValidatedFloat>;
template class Polynomial<UpperInterval>;

template Polynomial<ExactFloat>::Polynomial(SizeType);
template Void Polynomial<ExactFloat>::append(MultiIndex const&, ExactFloat const&);
template OutputStream& Polynomial<ExactFloat>::_write(OutputStream&) const;
template OutputStream& Polynomial<ExactFloat>::_write(OutputStream&, List<String> const&) const;

} //namespace Ariadne


