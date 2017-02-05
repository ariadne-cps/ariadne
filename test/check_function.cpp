/***************************************************************************
 *            check_function.cpp
 *
 *  Copyright 2009--17  Pieter Collins
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

#include "function/functional.hpp"

#include "function/function.hpp"
#include "function/taylor_model.hpp"
#include "geometry/box.hpp"
#include "algebra/algebra.hpp"

#include "numeric/float.decl.hpp"

#include "test.hpp"
#include "check_function.hpp"

using namespace Ariadne;

int main() {
    ARIADNE_CURRENT_TESTING_CLASS="EffectiveFunction";
    CheckFunctionConcept<EffectiveFunction>().check();
    return ARIADNE_TEST_FAILURES;
}

