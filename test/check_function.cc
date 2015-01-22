/***************************************************************************
 *            check_function.cc
 *
 *  Copyright 2009-14  Pieter Collins
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

#include "function/functional.h"

#include "function/function.h"
#include "geometry/box.h"

#include "numeric/float.decl.h"

#include "test.h"
#include "check_function.h"

using namespace Ariadne;

int main() {
    ARIADNE_CURRENT_TESTING_CLASS="EffectiveFunction";
    CheckFunctionConcept<EffectiveFunction>().check();
    return ARIADNE_TEST_FAILURES;
}

