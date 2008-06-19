/***************************************************************************
 *            numeric/integer.h
 *
 *  Copyright  2004-7  Alberto Casagrande, Pieter Collins
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
 
/*! \file numeric/integer.h
 *  \brief Multiple-precision integer type and interger functions.
 */

#ifndef ARIADNE_NUMERIC_INTEGER_H
#define ARIADNE_NUMERIC_INTEGER_H

//#include <gmpxx.h>
#include <gmp.h>
#include <iosfwd>
#include <stdexcept>

#include "numeric/traits.h"
#include "numeric/expression.h"
#include "numeric/operators.h"

#include "numeric/integer.class.h"
#include "numeric/integer.inline.h"
#include "numeric/integer.template.h"
 
namespace Ariadne {
  

    // Declare stream i/o operators
    std::ostream& operator<<(std::ostream& os, const Integer& n);
    std::istream& operator>>(std::istream& is, Integer& n);


} // namespace Ariadne

#endif /* ARIADNE_NUMERIC_INTEGER_H */
