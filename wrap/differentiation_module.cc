/***************************************************************************
 *            differentiation_module.cc
 *
 *  Copyright  2007-8  Pieter Collins
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

#include "boost_python.h"

namespace Ariadne {
class Float;
class ExactInterval;
template<class X> class DenseDifferential;
template<class X> class SparseDifferential;
template<class DIFF> class DifferentialVector;
}

template<class DIFF> Void export_differential();
template<class DIFF> Void export_differential_vector();

using namespace Ariadne;

BOOST_PYTHON_MODULE(differentiation)
{
    export_differential< DenseDifferential<Float> >();
    export_differential< DenseDifferential<ExactInterval> >();
    export_differential< SparseDifferential<Float> >();
    export_differential< SparseDifferential<ExactInterval> >();

    export_differential_vector< DenseDifferential<Float> >();
    export_differential_vector< DenseDifferential<ExactInterval> >();
    export_differential_vector< SparseDifferential<Float> >();
    export_differential_vector< SparseDifferential<ExactInterval> >();
}
